/**
 * Copyright 2022 Alban Fichet
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 *   3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <sstream>
#include <limits>
#include <regex>
#include <string>
#include <map>
#include <set>
#include <algorithm>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include <JXLImage.h>
#include <EXRSpectralImage.h>

#include <moments.h>
#include <moments_image.h>


void compress_spectral_framebuffer(
    const SpectralFramebuffer* framebuffer,
    std::vector<std::vector<float>>& compressed_moments,
    std::vector<float>& mins,
    std::vector<float>& maxs)
{
    std::vector<double> phases;
    std::vector<double> moments_image;
    std::vector<double> compressed_moments_image;
    
    const uint32_t n_moments = framebuffer->wavelengths_nm.size();
    const uint32_t n_pixels = framebuffer->image_data.size() / framebuffer->wavelengths_nm.size();

    // TODO: this shall be revisited
    std::vector<double> spectral_wavelengths(framebuffer->wavelengths_nm.size());
    std::vector<double> spectral_framebuffer(framebuffer->image_data.size());

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_wavelengths.size(); i++) {
        spectral_wavelengths[i] = framebuffer->wavelengths_nm[i];
    }

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_framebuffer.size(); i++) {
        double v = framebuffer->image_data[i];

        if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
            v = 1e-8;
        }

        spectral_framebuffer[i] = v;
    }

    wavelengths_to_phases(spectral_wavelengths, phases);
    
    compute_moments_image(
        phases, 
        spectral_framebuffer,
        n_pixels, 1, 
        n_moments, 
        moments_image
    );

    unbounded_to_bounded_compress_moments_image(
        moments_image,
        n_pixels, 1,
        n_moments,
        compressed_moments_image
    );

    compressed_moments.resize(n_moments);

    for (size_t m = 0; m < compressed_moments.size(); m++) {
        compressed_moments[m].resize(n_pixels);
    }

    // DC component does not need further modification
    for (size_t i = 0; i < n_pixels; i++) {
        compressed_moments[0][i] = compressed_moments_image[n_moments * i];
    }

    // Rescale AC components in [0..1]
    for (size_t m = 1; m < n_moments; m++) {
        // Get min / max
        double v_min = compressed_moments_image[m];
        double v_max = compressed_moments_image[m];

        for (size_t i = 0; i < n_pixels; i++) {
            v_min = std::min(v_min, compressed_moments_image[n_moments * i + m]);
            v_max = std::max(v_max, compressed_moments_image[n_moments * i + m]);
        }

        mins.push_back(v_min);
        maxs.push_back(v_max);

        // Now rescale moments
        for (size_t i = 0; i < n_pixels; i++) {
            const float v = compressed_moments_image[n_moments * i + m];

            compressed_moments[m][i] = (v - v_min) / (v_max - v_min);
        }
    }
}


void quantization_from_exr(PixelType type, size_t& n_bits, size_t& n_exponent_bits) {
    switch (type) {
        case PixelType::UINT:
            n_bits = 32;
            n_exponent_bits = 0;
            break;

        case PixelType::HALF:
            n_bits = 16;
            n_exponent_bits = 5;
            break;

        case PixelType::FLOAT:
            n_bits = 32;
            n_exponent_bits = 8;
            break;

        default:
            throw std::runtime_error("Unknown pixel type");
            break;
    }
}


int main(int argc, char *argv[]) 
{
    // TODO: generate a default output name and check if exists,
    // this will allow drag and drop of an EXR over the executable
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_in> <jxl_out>" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    EXRSpectralImage image_in(filename_in);

    const std::vector<SpectralFramebuffer*>& spectral_framebuffers = image_in.getSpectralFramebuffers();
    const std::vector<GreyFramebuffer*>& extra_framebuffers = image_in.getExtraFramebuffers();

    // get_buffers_from_exr(exr_in, spectral_framebuffers, extra_framebuffers);

    JXLImage jxl_out(image_in.width(), image_in.height());
    SGEGBox box;

    box.exr_attributes = image_in.getAttributesData();
    
    for (const SpectralFramebuffer* fb: spectral_framebuffers) {
        SGEGSpectralGroup sg;

        sg.root_name.resize(fb->root_name.size() + 1);
        std::memcpy(sg.root_name.data(), fb->root_name.c_str(), sg.root_name.size() * sizeof(char));
        sg.wavelengths = fb->wavelengths_nm;

        size_t main_n_bits;
        size_t main_n_exponent_bits;

        quantization_from_exr(fb->pixel_type, main_n_bits, main_n_exponent_bits);

        std::cout << main_n_bits << " " << main_n_exponent_bits << std::endl;

        std::vector<std::vector<float>> compressed_moments;

        compress_spectral_framebuffer(fb, compressed_moments, sg.mins, sg.maxs);

        // Now we can save to JPEG XL
        for (size_t m = 0; m < compressed_moments.size(); m++) {
            float n_bits;
            float n_exponent_bits;

            if (m == 0) {
                n_bits          = main_n_bits;
                n_exponent_bits = main_n_exponent_bits;
            } else {
                n_bits          = 8;
                n_exponent_bits = 0;
            }

            const size_t idx = jxl_out.appendFramebuffer(
                compressed_moments[m],
                1,
                n_bits,
                n_exponent_bits,
                1,
                fb->root_name.c_str());            
            
            sg.layer_indices.push_back(idx);
        }

        box.spectral_groups.push_back(sg);
    }

    for (const GreyFramebuffer* fb: extra_framebuffers) {
        SGEGGrayGroup gg;
        
        gg.layer_name.resize(fb->layer_name.size() + 1);
        std::memcpy(gg.layer_name.data(), fb->layer_name.c_str(), gg.layer_name.size() * sizeof(char));

        size_t n_bits;
        size_t n_exponent_bits;

        quantization_from_exr(fb->pixel_type, n_bits, n_exponent_bits);

        gg.layer_index = jxl_out.appendFramebuffer(
            fb->image_data,
            1,
            n_bits,
            n_exponent_bits,
            1,
            fb->layer_name.c_str()
        );

        box.gray_groups.push_back(gg);
    }

    // box.print();

    jxl_out.setBox(box);
    jxl_out.write(filename_out);

    return 0;
}
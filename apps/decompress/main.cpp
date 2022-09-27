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
#include <cassert>
#include <algorithm>

#include <JXLImage.h>
#include <EXRSpectralImage.h>
// TODO: remove
// #include <EXRImage.h>

#include <moments_image.h>
#include <moments.h>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>


void decompress_spectral_framebuffer(
    const std::vector<float>& wavelengths,
    const std::vector<std::vector<float>>& compressed_moments,
    const std::vector<float>& mins,
    const std::vector<float>& maxs,
    std::vector<float>& spectral_framebuffer)
{
    // const uint32_t n_pixels = compressed_moments.size() / wavelengths.size();
    std::vector<float> phases;
    std::vector<float> moments_image;

    const size_t n_moments = compressed_moments.size();
    const size_t n_pixels = compressed_moments[0].size();
    
    wavelengths_to_phases(wavelengths, phases);

    std::vector<float> compressed_moments_rescaled(n_moments * n_pixels);

    for (size_t i = 0; i < n_pixels; i++) {
        compressed_moments_rescaled[n_moments * i] = compressed_moments[0][i];
    }

    for (size_t m = 1; m < n_moments; m++) {
        const float v_min = mins[m - 1];
        const float v_max = maxs[m - 1];

        for (size_t i = 0; i < n_pixels; i++) {
            const float v = compressed_moments[m][i];

            compressed_moments_rescaled[n_moments * i + m] = (v_max - v_min) * v + v_min;
        }
    }

    unbounded_decompress_moments_image(compressed_moments_rescaled, n_pixels, 1, n_moments, moments_image);
    compute_density_image(phases, moments_image, n_pixels, 1, n_moments, spectral_framebuffer);
}


PixelType quantization_to_exr(size_t n_bits, size_t n_exponent_bits) {
    if (n_exponent_bits > 0) {
        if (n_bits > 16) {
            return PixelType::FLOAT;
        } else {
            return PixelType::HALF;
        }
    } else {
        return PixelType::UINT;
    }
}



int main(int argc, char* argv[]) 
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <jxl_in> <exr_out>" << std::endl;
        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    const JXLImage jxl_image(filename_in);

    const SGEGBox box   = jxl_image.getBox();
    const size_t width  = jxl_image.width();
    const size_t height = jxl_image.height();

    EXRSpectralImage exr_out(width, height);

    exr_out.setAttributesData(box.exr_attributes);

    for (const SGEGSpectralGroup& sg: box.spectral_groups) {
        std::string root_name = sg.root_name.data();
        const size_t n_moments = sg.layer_indices.size();

        assert(sg.layer_indices.size() == (sg.mins.size() + 1));
        assert(sg.mins.size() == sg.maxs.size());

        std::vector<std::vector<float>> moments(n_moments);
        std::vector<float> spectral_framebuffer;

        for (size_t m = 0; m < n_moments; m++) {
            moments[m] = jxl_image.getFramebufferDataConst(sg.layer_indices[m]);
        }

        decompress_spectral_framebuffer(
            sg.wavelengths, 
            moments, 
            sg.mins, sg.maxs, 
            spectral_framebuffer
        );
        
        const JXLFramebuffer* main_fb = jxl_image.getFramebuffer(sg.layer_indices[0]);
        const size_t n_bits           = main_fb->getBitsPerSample();
        const size_t n_exponent_bits  = main_fb->getExponentBitsPerSample();
        const PixelType pixel_type    = quantization_to_exr(n_bits, n_exponent_bits);

        exr_out.appendSpectralFramebuffer(
            sg.wavelengths,
            spectral_framebuffer,
            root_name,
            pixel_type
        );
    }

    for (const SGEGGrayGroup& gg: box.gray_groups) {
        const JXLFramebuffer* main_fb = jxl_image.getFramebuffer(gg.layer_index);
        const size_t n_bits           = main_fb->getBitsPerSample();
        const size_t n_exponent_bits  = main_fb->getExponentBitsPerSample();
        const PixelType pixel_type    = quantization_to_exr(n_bits, n_exponent_bits);

        exr_out.appendExtraFramebuffer(
            jxl_image.getFramebufferDataConst(gg.layer_index), 
            gg.layer_name.data(),
            pixel_type
        );
    }

    exr_out.write(filename_out);

    return 0;
}

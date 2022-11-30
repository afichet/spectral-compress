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
    const std::vector<uint8_t>& relative_scales,
    float global_max,
    std::vector<float>& spectral_framebuffer)
{
    const size_t n_moments = compressed_moments.size();
    const size_t n_pixels = compressed_moments[0].size();

    std::vector<double> wavelengths_d(wavelengths.size());
    std::vector<double> compressed_moments_d(n_moments * n_pixels);
    std::vector<double> mins_d(n_moments - 1), maxs_d(n_moments - 1);
    std::vector<double> spectral_framebuffer_d;

    for (size_t i = 0; i < wavelengths.size(); i++) {
        wavelengths_d[i] = wavelengths[i];
    }

    for (size_t px = 0; px < n_pixels; px++) {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments_d[px * n_moments + m] = compressed_moments[m][px];
        }
    }

    for (size_t i = 0; i < n_moments - 1; i++) {
        mins_d[i] = mins[i];
        maxs_d[i] = maxs[i];
    }

    upperbound_decompress_spectral_image(
        wavelengths_d,
        compressed_moments_d,
        mins_d, maxs_d,
        relative_scales,
        global_max,
        n_pixels, n_moments,
        spectral_framebuffer_d
    );

    spectral_framebuffer.resize(spectral_framebuffer_d.size());

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_framebuffer_d.size(); i++) {
        spectral_framebuffer[i] = spectral_framebuffer_d[i];
    }
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
        // 1 extra framebuffer is used to store relative max
        const size_t n_moments = sg.layer_indices.size() - 1;

        assert(sg.layer_indices.size() == (sg.mins.size() + 1));
        assert(sg.mins.size() == sg.maxs.size());

        std::vector<std::vector<float>> moments(n_moments);
        std::vector<float> relative_scales_f;
        std::vector<float> spectral_framebuffer;

        for (size_t m = 0; m < n_moments; m++) {
            moments[m] = jxl_image.getFramebufferDataConst(sg.layer_indices[m]);
        }

        relative_scales_f = jxl_image.getFramebufferDataConst(sg.layer_indices.back());

        // TODO: This is a bit hacky
        std::vector<uint8_t> relative_scales(relative_scales_f.size());

        for (size_t i = 0; i < relative_scales_f.size(); i++) {
            relative_scales[i] = std::numeric_limits<uint8_t>::max() * relative_scales_f[i];
        }

        decompress_spectral_framebuffer(
            sg.wavelengths,
            moments,
            sg.mins, sg.maxs,
            relative_scales,
            sg.global_max,
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

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

#include <EXRSpectralImage.h>
#include <JXLImage.h>
#include <SGEG_box.h>

#include <moments.h>
#include <moments_image.h>


int main(int argc, char *argv[]) 
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_in> <jxl_out>" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    const SEXR::EXRSpectralImage image_in(filename_in);

    const size_t width     = image_in.width();
    const size_t height    = image_in.height();
    const size_t n_bands   = image_in.nSpectralBands();
    const size_t n_moments = n_bands;

    std::vector<float> phases;
    std::vector<float> wavelengths;
    
    for (size_t i = 0; i < n_bands; i++) {
        wavelengths.push_back(image_in.wavelength_nm(i));
    }

    wavelengths_to_phases(wavelengths, phases);

    SGEG_box sgeg_box(n_moments);
    sgeg_box.n_wl_original = n_bands;
    sgeg_box.wl_min_nm = wavelengths.front();
    sgeg_box.wl_max_nm = wavelengths.back();

    if (image_in.isEmissive()) {
        std::vector<float> moments_image(width * height * (n_moments + 1));
        std::vector<float> compressed_moments_image(width * height * (n_moments + 1));
        
        const float* framebuffer = &image_in.emissive(0, 0, 0, 0);
        
        compute_moments_image(
            phases.data(), phases.size(),
            framebuffer,
            width, height,
            n_moments,
            moments_image.data()
        );

        compress_moments_image(
            moments_image,
            width, height,
            n_moments,
            compressed_moments_image
        );

        std::vector<float> dc_component(width * height);
        // std::vector<std::vector<float>> rescaled_ac(n_moments);
        std::vector<std::vector<uint8_t>> rescaled_ac(n_moments);

        // DC component
        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            dc_component[i] = compressed_moments_image[(n_moments + 1) * i];
        }


        // AC components
        for (size_t m = 0; m < n_moments; m++) {
            // We need min and max of a given moment to later rescale in 0..1 for
            // integer quantization
            float min = 1999, max = 0;

            for (size_t i = 0; i < width * height; i++) {
                const float v = compressed_moments_image[(n_moments + 1) * i + m + 1];
                min = std::min(min, v);
                max = std::max(max, v);
            }

            sgeg_box.moment_min[m] = min;
            sgeg_box.moment_max[m] = max;

            rescaled_ac[m].resize(width * height);

            #pragma omp parallel for
            for (size_t i = 0; i < width * height; i++) {
                const float v = (compressed_moments_image[(n_moments + 1) * i + m + 1] - min) / (max - min);
                //rescaled_ac[m][i] = v;
                rescaled_ac[m][i] = 255.f * v;
            }
        }

        // Create the JXL
        JXLImageWriter jxl_image(width, height, sgeg_box, n_moments);

        jxl_image.addMainFramebuffer(dc_component.data());

        for (size_t m = 0; m < n_moments; m++) {
            jxl_image.addSubFramebuffer(rescaled_ac[m].data(), m);
        }

        jxl_image.save(filename_out);
    } else {
        std::cerr << "Only emissive images are supported now." << std::endl;
    }

    if (image_in.isReflective()) {
        std::cout << "[WARN] Reflective layers will be ignored." << std::endl;
    }

    if (image_in.isBispectral()) {
        std::cout << "[WARN] Bispectral layers will be ignored." << std::endl;
    }
    
    if (image_in.isPolarised()) {
        std::cout << "[WARN] The polarised layers will be ignored." << std::endl;
    }

    return 0;
}
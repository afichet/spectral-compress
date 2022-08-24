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

#include <limits>

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

    const bool saveRGB = false;

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    const SEXR::EXRSpectralImage image_in(filename_in);

    std::vector<float> rgb_image;
    image_in.getRGBImage(rgb_image);

    const size_t width     = image_in.width();
    const size_t height    = image_in.height();
    const size_t n_bands   = image_in.nSpectralBands();
    const size_t n_moments = n_bands - 1;

    std::vector<float> phases;
    std::vector<float> wavelengths;
    
    for (size_t i = 0; i < n_bands; i++) {
        wavelengths.push_back(image_in.wavelength_nm(i));
    }

    wavelengths_to_phases(wavelengths, phases);

    SGEG_box sgeg_box(n_moments, n_bands);
    sgeg_box.wavelengths = wavelengths;

    if (image_in.isEmissive()) {
        sgeg_box.is_reflective = false;

        std::vector<float> moments_image(width * height * (n_moments + 1));
        std::vector<float> compressed_moments_image(width * height * (n_moments + 1));
        
        const float* og_framebuffer = &image_in.emissive(0, 0, 0, 0);
        std::vector<float> framebuffer(width * height * n_bands);

        memcpy(framebuffer.data(), og_framebuffer, width * height * n_bands * sizeof(float));


        std::cout << "Conditioning framebuffers...\t";
        // Ensure positive definite values
        // We need to ensure as well that we do not have any spectrum with all
        // null values, otherwise this would break later at the compression of
        // moments.
        for (size_t i = 0; i < width * height * n_bands; i++) {
            const float v = framebuffer[i];

            // if (std::isinf(v) || std::isnan(v) || v <= 0) {
            //     framebuffer[i] = std::numeric_limits<float>::epsilon();
            // }


            if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
                framebuffer[i] = 1e-8;
            }
        }

        std::cout << "[DONE]" << std::endl;


        std::cout << "Computing moments...\t\t";

        compute_moments_image(
            phases,
            framebuffer,
            width, height,
            n_moments,
            moments_image
        );

        for (size_t i = 0; i < width * height; i++) {
            for (size_t m = 0; m < n_moments + 1; m++) {
                const float v = moments_image[(n_moments + 1) * i + m];
                if (std::isinf(v) || std::isnan(v)) {
                    std::cout << "Invalid moment value!" << std::endl;
                }
            }
        }

        compress_moments_image(
            moments_image,
            width, height,
            n_moments,
            compressed_moments_image
        );


        for (size_t i = 0; i < width * height; i++) {
            for (size_t m = 0; m < n_moments + 1; m++) {
                const float v = compressed_moments_image[(n_moments + 1) * i + m];
                if (std::isinf(v) || std::isnan(v)) {
                    std::cout << "Invalid compressed moment value!" << std::endl;
                }
            }
        }

        std::cout << "[DONE]" << std::endl;

        std::cout << "Rescaling moments...\t\t";


        std::vector<float> dc_component(width * height);
        std::vector<std::vector<float>> rescaled_ac(n_moments);

        // DC component: it bascially contains all the HDR information so
        // no further transformation is performed on this component.
        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            dc_component[i] = compressed_moments_image[(n_moments + 1) * i];
        }

        // AC components: those components are higher frequency and as such
        // can be more aggressively compressed. We reacale each component
        // to make use of the maximum range of an uint8.
        for (size_t m = 0; m < n_moments; m++) {
            // We need min and max of a given moment to later rescale in 0..1 for
            // integer quantization
            float min = compressed_moments_image[(n_moments + 1) + m + 1];
            float max = compressed_moments_image[(n_moments + 1) + m + 1];

            for (size_t i = 0; i < width * height; i++) {
                const float v = compressed_moments_image[(n_moments + 1) * i + m + 1];
                if (!std::isinf(v)) {
                    min = std::min(min, v);
                    max = std::max(max, v);
                }
            }

            // TODO: FIX THE CLAMPING
            sgeg_box.moment_min[m] = std::max(-1.f, min);
            sgeg_box.moment_max[m] = std::min(1.f, max);

            rescaled_ac[m].resize(width * height);

            #pragma omp parallel for
            for (size_t i = 0; i < width * height; i++) {
                const float v = (compressed_moments_image[(n_moments + 1) * i + m + 1] - min) / (max - min);

                // TODO: FIX THE CLAMPING
                rescaled_ac[m][i] = std::max(0.f, std::min(1.f, v));
            }
        }

        std::cout << "[DONE]" << std::endl;

        std::cout << "Adding framebuffers to JXL...\t";

        // Debug
        sgeg_box.print_info();

        // Create the JXL
        JXLImageWriter jxl_image(width, height, sgeg_box, saveRGB ? n_moments + 1 : n_moments);

        if (saveRGB) {
            // Convert to RGB 24 bits
            std::vector<uint8_t> rgb_image_ldr(width * height * 3);

            for (size_t i = 0; i < width * height; i++) {
                for (size_t c = 0; c < 3; c++) {
                    rgb_image_ldr[3 * i + c] = 255.f * std::min(1.f, std::max(0.f, rgb_image[3 * i + c]));
                }
            }

            jxl_image.addMainRGBFramebuffer(rgb_image_ldr.data());

            jxl_image.addSubFramebuffer(dc_component.data(), 0, 1, true);

            for (size_t m = 0; m < n_moments; m++) {
                uint32_t bits_per_pixel = 8;

                if (m >= 4) bits_per_pixel = 6;
                if (m >= 6) bits_per_pixel = 4;

                jxl_image.addSubFramebuffer(rescaled_ac[m].data(), m + 1, 1, false, bits_per_pixel);
            }
        } else {
            jxl_image.addMainFramebuffer(dc_component.data());

            for (size_t m = 0; m < n_moments; m++) {
                uint32_t bits_per_pixel = 8;

                if (m >= 4) bits_per_pixel = 6;
                if (m >= 6) bits_per_pixel = 4;

                jxl_image.addSubFramebuffer(rescaled_ac[m].data(), m, 1, false, bits_per_pixel);
            }
        }

        std::cout << "[DONE]" << std::endl;

        std::cout << "Saving file...\t\t\t";

        jxl_image.save(filename_out);

        std::cout << "[DONE]" << std::endl;

    } else {
        sgeg_box.is_reflective = true;

        std::vector<float> reflective_fb(width * height * n_bands);
        std::vector<float> moments_image(width * height * (n_moments + 1));
        std::vector<float> compressed_bounded_moments_image(width * height * (n_moments + 1));

        memcpy(reflective_fb.data(), &image_in.reflective(0, 0, 0), reflective_fb.size() * sizeof(float));

        for (size_t i = 0; i < width * height * n_bands; i++) {
            const float v = reflective_fb[i];

            if (std::isinf(v) || std::isnan(v) || v < 0) {
                reflective_fb[i] = 0;
            }
        }

        compute_moments_image(
            phases,
            reflective_fb,
            width, height,
            n_moments,
            moments_image
        );

        compress_bounded_moments_image(
            moments_image,
            width, height,
            n_moments,
            compressed_bounded_moments_image
        );

        std::vector<float> dc_component(width * height);
        // std::vector<std::vector<uint8_t>> rescaled_ac(n_moments);
        std::vector<std::vector<float>> rescaled_ac(n_moments);

        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            dc_component[i] = compressed_bounded_moments_image[(n_moments + 1) * i];
        }


        for (size_t m = 0; m < n_moments; m++) {
            float min = 1.;
            float max = -1.;


            for (size_t i = 0; i < width * height; i++) {
                const float v = compressed_bounded_moments_image[(n_moments + 1) * i + m + 1];
                if (!std::isinf(v)) {
                    min = std::min(min, v);
                    max = std::max(max, v);
                }
            }

            sgeg_box.moment_min[m] = min;
            sgeg_box.moment_max[m] = max;

            rescaled_ac[m].resize(width * height);

            #pragma omp parallel for
            for (size_t i = 0; i < width * height; i++) {
                const float v = (compressed_bounded_moments_image[(n_moments + 1) * i + m + 1] - min) / (max - min);
                // rescaled_ac[m][i] = 255.f * v;
                rescaled_ac[m][i] = v;
            }
        }

        JXLImageWriter jxl_image(width, height, sgeg_box, n_moments);
        jxl_image.addMainFramebuffer(dc_component.data());

        for (size_t m = 0; m < n_moments; m++) {
            jxl_image.addSubFramebuffer(rescaled_ac[m].data(), m);
        }

        jxl_image.save(filename_out);
    }


    if (image_in.isBispectral()) {
        std::cout << "[WARN] Bispectral layers will be ignored." << std::endl;
    }
    
    if (image_in.isPolarised()) {
        std::cout << "[WARN] The polarised layers will be ignored." << std::endl;
    }

    return 0;
}
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

#include <EXRSpectralImage.h>
#include <JXLImage.h>
#include <SGEG_box.h>

#include <moments_image.h>
#include <moments.h>


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

    const JXLImageReader jxl_image(filename_in);
    const SGEG_box sgeg_box  = jxl_image.get_sgeg();
    const size_t width       = jxl_image.width();
    const size_t height      = jxl_image.height();
    const uint32_t n_moments = sgeg_box.n_moments;

    // Debug
    // jxl_image.print_basic_info();
    // sgeg_box.print_info();

    std::vector<float> dc_image(width * height);
    std::vector<std::vector<float>> ac_images(sgeg_box.n_moments);

    jxl_image.getMainFramebuffer(dc_image);

    for (uint32_t m = 0; m < n_moments; m++) {
        jxl_image.getSubFramebuffer(ac_images[m], m);

        // Apply rescaling
        for (size_t i = 0; i < width * height; i++) {
            ac_images[m][i] = ac_images[m][i] * (sgeg_box.moment_max[m] - sgeg_box.moment_min[m]) + sgeg_box.moment_min[m];
        }
    }

    std::vector<float> compressed_moments_image(width * height * (n_moments + 1));
    std::vector<float> moments_image(width * height * (n_moments + 1));

    for (size_t i = 0; i < width * height; i++) {
        compressed_moments_image[(n_moments + 1) * i] = dc_image[i]; 
    }

    for (uint32_t m = 0; m < n_moments; m++) {
        for (size_t i = 0; i < width * height; i++) {
            compressed_moments_image[(n_moments + 1) * i + m + 1] = ac_images[m][i]; 
        }
    }

    decompress_moments_image(
        compressed_moments_image,
        width, height, n_moments,
        moments_image
    );

    const std::vector<float> wavelengths_nm = sgeg_box.wavelengths;
    std::vector<float> phases;

    wavelengths_to_phases(wavelengths_nm, phases);

    SEXR::EXRSpectralImage image_out(width, height, wavelengths_nm, SEXR::SpectrumType::EMISSIVE);

    compute_density_image(
        phases.data(), phases.size(),
        moments_image.data(),
        width,
        height,
        n_moments,
        &image_out.emissive(0, 0, 0, 0)
    );

    image_out.save(filename_out);

    return 0;
}

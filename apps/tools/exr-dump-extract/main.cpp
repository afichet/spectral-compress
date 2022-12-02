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
#include <cmath>
#include <vector>
#include <lodepng.h>
#include <EXRSpectralImage.h>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <exr-dump>" << std::endl;

        return 0;
    }

    const float exposure = 5.f;
    const float gamma = 1.f/2.2f;

    EXRSpectralImage* img = EXRSpectralImage::read_dump(argv[1]);

    for(const SpectralFramebuffer* fb: img->getSpectralFramebuffers()) {
        for (size_t i = 0; i < fb->wavelengths_nm.size(); i++) {
            std::stringstream filename;
            filename << fb->root_name << "." << fb->wavelengths_nm[i] << ".png";

            std::vector<uint8_t> out_fb(4 * img->width() * img->height());

            for (size_t p = 0; p < img->width() * img->height(); p++) {
                const float pixel = fb->image_data[p * fb->wavelengths_nm.size() + i];
                const float value = std::max(0.f, std::min(1.f, std::pow(std::exp2(exposure) * pixel, gamma)));

                out_fb[4 * p + 0] = 255 * value;
                out_fb[4 * p + 1] = 255 * value;
                out_fb[4 * p + 2] = 255 * value;
                out_fb[4 * p + 3] = 255;
            }

            lodepng::encode(filename.str(), out_fb.data(), img->width(), img->height());
        }
    }

    delete img;

    return 0;
}

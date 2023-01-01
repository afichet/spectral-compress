/**
 * Copyright 2022 - 2023 Alban Fichet
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

#include <EXRSpectralImage.h>
#include "macbeth_data.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>


int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <EXR image out>" << std::endl;

        exit(0);
    }

    const size_t patch_sz = 50;

    const size_t n_patches_x = 6;
    const size_t n_patches_y = 4;

    const size_t width  = n_patches_x * patch_sz;
    const size_t height = n_patches_y * patch_sz;

    const size_t n_bands = 36;

    EXRSpectralImage image_out(width, height);

    std::vector<float> wavelengths(n_bands);
    std::memcpy(wavelengths.data(), macbeth_wavelengths, n_bands * sizeof(float));

    std::vector<float> framebuffer(width * height * n_bands);

    for (size_t y = 0; y < height; y++) {
        const size_t v = (y * n_patches_y) / height;
        assert(v < n_patches_y);

        for (size_t x = 0; x < width; x++) {
            const size_t u = (x * n_patches_x) / width;
            assert(u < n_patches_x);

            std::memcpy(
                &framebuffer[n_bands * (y * width + x)],
                macbeth_patches[v * n_patches_x + u],
                n_bands * sizeof(float)
            );
        }
    }

    image_out.appendSpectralFramebuffer(wavelengths, framebuffer, "T");

    image_out.write(argv[1]);

    return 0;
}

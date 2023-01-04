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

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include <JXLImage.h>

#include <sstream>

/**
 * Performs various tests on the JXL compressor
 * - Quantization: we expect progressively high banding
 * - Downsampling: we expect increasingly lower resolution
 */
int main(int argc, char* argv[])
{
    (void)argc; (void)argv;

    size_t width = 640;
    size_t height = 480;

    std::vector<float> framebuffer(width * height);

    // ------------------------------------------------------------------------
    // Generate a gradiant to check if the quantization ratio works as expected
    // ------------------------------------------------------------------------

    for (size_t y = 0; y < height; y++) {
        const float value = (float)y / (float)(height - 1);

        for (size_t x = 0; x < width; x++) {
            framebuffer[y * width + x] = value;
        }
    }

    // Save multiple files with a different quantization
    for (size_t q = 8; q > 1; q--) {
        std::stringstream ss;
        ss << "quantization_" << q << ".jxl";

        JXLImage jxl_out(width, height);
        jxl_out.appendFramebuffer(
            framebuffer, 1,
            q, 0,
            1, 0.f
        );

        jxl_out.write(ss.str());
    }

    // ------------------------------------------------------------------------
    // Generate lines of different heigth to check if the downsampling works
    // as expected
    // ------------------------------------------------------------------------

    int start_sz = 1;
    int remaining_lines = start_sz;
    int curr_sz = start_sz;

    for (size_t y = 0; y < height; y++) {
        float value = 0;

        if (remaining_lines > 0) {
            value = 1;
        } else if (-remaining_lines < curr_sz) {
            value = 0;
        } else {
            curr_sz += 1;
            remaining_lines = curr_sz;
        }

        for (size_t x = 0; x < width; x++) {
            framebuffer[y * width + x] = value;
        }

        --remaining_lines;
    }

    // It seems JXL supports this downsampling rates only
    size_t downsampling_ratios[] = {1, 2, 4, 8};

    for (size_t downsampling: downsampling_ratios) {
        std::stringstream ss;
        ss << "downsampled_" << downsampling << ".jxl";

        JXLImage jxl_out(width, height);
        jxl_out.appendFramebuffer(
            framebuffer, 1,
            8, 0,
            downsampling, 0.f
        );

        jxl_out.write(ss.str());
    }

    return 0;
}

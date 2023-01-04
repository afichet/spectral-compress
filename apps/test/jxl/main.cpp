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

#include "lodepng.h"


void write_foo_image(const char* filename, size_t width, size_t height)
{
    // Generate an image
    std::vector<float> main_framebuffer(3 * width * height);
    std::vector<float> sub_framebuffer(width * height);

    for (size_t y = 0; y < height; y++) {
        const float v = float(y) / float(height - 1);

        for (size_t x = 0; x < width; x++) {
            const float u = float(x) / float(width - 1);

            main_framebuffer[3 * (y * width + x) + 0] = u;
            main_framebuffer[3 * (y * width + x) + 1] = u*v;
            main_framebuffer[3 * (y * width + x) + 2] = 0;

            sub_framebuffer [y * width + x] = v;
        }
    }

    // Output to JXL
    JXLImage image_out(width, height);

    image_out.appendFramebuffer(main_framebuffer, 3);
    image_out.appendFramebuffer(sub_framebuffer, 1, 8, 0, 1, "Bonjour!");

    image_out.write(filename);
}


void to_sRGBA(
    const std::vector<float>& in_framebuffer,
    uint32_t n_color_channels,
    std::vector<uint8_t>& out_framebuffer)
{
    out_framebuffer.resize(in_framebuffer.size() * 4);

    switch (n_color_channels) {
        case 1:
            for (size_t i = 0; i < in_framebuffer.size(); i++) {
                for (int c = 0; c < 3; c++) {
                    out_framebuffer[4 * i + c] = 255 * std::pow(in_framebuffer[i], 1.f/2.2f);
                }

                out_framebuffer[4 * i + 3] = 255;
            }
            break;
        case 3:
            for (size_t i = 0; i < in_framebuffer.size(); i++) {
                for (int c = 0; c < 3; c++) {
                    out_framebuffer[4 * i + c] = 255 * std::pow(in_framebuffer[3 * i + c], 1.f/2.2f);
                }

                out_framebuffer[4 * i + 3] = 255;
            }
            break;
        default:
            break;
    }
}


int main(int argc, char* argv[])
{
    (void)argc; (void)argv;

    write_foo_image("Hello.jxl", 640, 480);

    JXLImage jxl_image("Hello.jxl");

    const uint32_t width = jxl_image.width();
    const uint32_t height = jxl_image.height();

    for (size_t i = 0; i < jxl_image.n_framebuffers(); i++) {
        std::stringstream ss;
        ss << "sub_" << i << ".png";

        const std::vector<float>& fb = jxl_image.getFramebufferData(i);
        const uint32_t n_color_channels = jxl_image.getFramebuffer(i)->getNColorChannels();

        std::vector<uint8_t> png;
        to_sRGBA(fb, n_color_channels, png);

        lodepng::encode(ss.str().c_str(), png.data(), width, height);
    }

    return 0;
}

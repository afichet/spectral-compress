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
#include <vector>
#include <sstream>
#include <cmath>

#include <JXLImage.h>


#include <jxl/encode.h>
#include <jxl/encode_cxx.h>

#include <jxl/thread_parallel_runner.h>
#include <jxl/thread_parallel_runner_cxx.h>

#include "lodepng.h"


void write_foo_image(const char* filename, size_t width, size_t height)
{
    // Generate an image
    std::vector<float> main_framebuffer(width * height);
    std::vector<uint8_t> sub_framebuffer(width * height);

    for (size_t y = 0; y < height; y++) {
        const float v = float(y) / float(height - 1);

        for (size_t x = 0; x < width; x++) {
            const float u = float(x) / float(width - 1);

            main_framebuffer[y * width + x] = u;
            sub_framebuffer [y * width + x] = 255 * v;
        }
    }

    SGEG_box sgeg(2);

    // Output to JXL
    JXLImageWriter w(width, height, sgeg, 1);

    w.addMainFramebuffer(main_framebuffer.data());
    w.addSubFramebuffer(sub_framebuffer.data(), 0, 1, 8, "Bonjour!");

    w.save(filename);
}


void to_sRGBA(
    const std::vector<float>& in_framebuffer,
    std::vector<uint8_t>& out_framebuffer)
{
    out_framebuffer.resize(in_framebuffer.size() * 4);

    for (size_t i = 0; i < in_framebuffer.size(); i++) {
        for (int c = 0; c < 3; c++) {
            out_framebuffer[4 * i + c] = 255 * std::pow(in_framebuffer[i], 1.f/2.2f);
        }

        out_framebuffer[4 * i + 3] = 255;
    }
}


int main(int argc, char* argv[])
{
    (void)argc; (void)argv;
    
    // write_foo_image("toto.jxl", 640, 480);

    // JXLImageReader r("toto.jxl");
    JXLImageReader r("Hello.jxl");

    size_t r_width = r.width();
    size_t r_height = r.height();

    r.print_basic_info();

    std::vector<float> main_image;
    r.getMainFramebuffer(main_image);

    std::vector<uint8_t> png_image;
    to_sRGBA(main_image, png_image);

    lodepng::encode("main.png", png_image.data(), r_width, r_height);

    for (size_t i = 0; i < r.n_subframebuffers(); i++) {
        std::stringstream ss;
        ss << "sub_" << i << ".png";

        std::vector<float> fb;
        r.getSubFramebuffer(fb, i);

        std::vector<uint8_t> png;
        to_sRGBA(fb, png);

        lodepng::encode(ss.str().c_str(), png.data(), r_width, r_height);
    }

    return 0;
}
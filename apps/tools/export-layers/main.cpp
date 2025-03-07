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
#include <sstream>
#include <cstdint>
#include <cstring>
#include <vector>

#include <lodepng.h>

#include <JXLImage.h>
#include <Util.h>

void export_float_gray_to_png(
    const std::vector<float>& framebuffer,
    uint32_t width, uint32_t height,
    const char* filename)
{
    std::vector<uint8_t> export_fb(framebuffer.size() * 4);

    for (size_t i = 0; i < framebuffer.size(); i++) {
        const uint8_t value = 255.f * Util::clamp(framebuffer[i], 0.f, 1.f);

        for (size_t c = 0; c < 3; c++) {
            export_fb[4 * i + c] = value;
        }

        export_fb[4 * i + 3] = 255;
    }

    lodepng::encode(filename, export_fb.data(), width, height);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <jxl_in>" << std::endl;
        exit(0);
    }

    JXLImage image_in(argv[1]);

    const uint32_t width  = image_in.width();
    const uint32_t height = image_in.height();

    // Export the main framebuffer
    std::vector<float> main_fb(width * height);
    std::memcpy(main_fb.data(), image_in.getFramebufferData(0).data(), width * height * sizeof(float));

    float min_v = main_fb[0], max_v = main_fb[0];

    // Rescaling
    for (uint32_t i = 0; i < main_fb.size(); i++) {
        min_v = std::min(min_v, main_fb[i]);
        max_v = std::max(max_v, main_fb[i]);
    }

    for (uint32_t i = 0; i < main_fb.size(); i++) {
        main_fb[i] = 1000 * main_fb[i] * (max_v - min_v) + min_v;
    }

    std::stringstream ss;
    ss << argv[1] << "_0.png";

    export_float_gray_to_png(main_fb, width, height, ss.str().c_str());


    for (uint32_t layer = 1; layer < image_in.n_framebuffers(); layer++) {
        std::stringstream ss;
        ss << argv[1] << "_" << layer << ".png";

        export_float_gray_to_png(
            image_in.getFramebufferData(layer),
            width, height,
            ss.str().c_str());
    }

    return 0;
}

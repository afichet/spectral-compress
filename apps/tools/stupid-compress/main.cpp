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

#include <cstring>

#include <EXRImage.h>
#include <JXLImage.h>

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl;
        std::cout << "------" << std::endl;
        std::cout << argv[0] << " <exr_in> <jxl_out> [<framedistance> = 0.1]" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];
    const float framedistance = (argc >= 4) ? std::stof(argv[3]) : 0.1f;

    const EXRImage exr_in(filename_in);
    JXLImage jxl_out(exr_in.width(), exr_in.height());
    SGEGBox box;

    box.exr_attributes = exr_in.getAttributesData();

    const std::vector<EXRFramebuffer*>& framebuffers = exr_in.getFramebuffersConst();

    for (const EXRFramebuffer* fb: framebuffers) {
        SGEGGrayGroup gg;
        const char* layer_name = fb->getName();

        gg.layer_index = jxl_out.appendFramebuffer(fb->getPixelDataConst(), 1, std::make_pair(32, 8), 1, framedistance, fb->getName());

        gg.layer_name.resize(std::strlen(layer_name) + 1);
        std::memcpy(gg.layer_name.data(), layer_name, gg.layer_name.size() * sizeof(char));

        box.gray_groups.push_back(gg);
    }

    jxl_out.setBox(box);

    std::cout << "Len ATTR: " << box.exr_attributes.size() << std::endl;

    jxl_out.write(filename_out);

    return 0;
}

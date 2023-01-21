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
#include <exception>

#include <EXRImage.h>
#include <JXLImage.h>

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl;
        std::cout << "------" << std::endl;
        std::cout << argv[0] << "<jxl_in> <exr_out>" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    try {
        const JXLImage jxl_in(filename_in);
        const SGEGBox& box = jxl_in.getBox();

        EXRImage exr_out(jxl_in.width(), jxl_in.height());

        for (size_t i = 0; i < box.gray_groups.size(); i++) {
            const SGEGGrayGroup& gg = box.gray_groups[i];

            exr_out.appendFramebuffer(
                jxl_in.getFramebufferDataConst(gg.layer_index),
                gg.layer_name.data()
            );
        }

        exr_out.setAttributesData(box.exr_attributes);

        exr_out.write(filename_out);
    } catch (std::exception &e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

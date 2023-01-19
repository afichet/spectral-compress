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
#include <string>
#include <cstring>

#include <tclap/CmdLine.h>

#include <EXRImage.h>

int main (int argc, char* argv[])
{
    std::string path_image_in, path_image_out;

    try {
        TCLAP::CmdLine cmd("Removes RGB layers from an OpenEXR image");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Path to OpenEXR image to remove RGB layers from.", true, "in.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Path to the OpenEXR image to write into.", true, "out.exr", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        cmd.parse(argc, argv);

        path_image_in = inputFileArg.getValue();
        path_image_out = outputFileArg.getValue();
    } catch (TCLAP::ArgException  &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        return 1;
    }

    const EXRImage image_in(path_image_in);

    EXRImage image_out(image_in.width(), image_in.height());

    for (const EXRFramebuffer* fb: image_in.getFramebuffersConst()) {
        if (std::strcmp(fb->getName(), "R") != 0
         && std::strcmp(fb->getName(), "G") != 0
         && std::strcmp(fb->getName(), "B") != 0) {
            image_out.appendFramebuffer(fb->getPixelDataConst(), fb->getName());
        }
    }

    image_out.setAttributesData(image_in.getAttributesData());

    image_out.write(path_image_out);
}

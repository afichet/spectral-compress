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
#include <exception>
#include <fstream>

#include <JXLImage.h>
#include <Util.h>

#include <tclap/CmdLine.h>


int main(int argc, char* argv[])
{
    std::string filename_in, filename_out;

    try {
        TCLAP::CmdLine cmd("Utility to compress all OpenEXR framebuffer in JPEG-XL.");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Specifies the OpenEXR input image.", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the JPEG-XL output image.", true, "output.jxl", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        cmd.parse(argc, argv);

        filename_in  = inputFileArg.getValue();
        filename_out = outputFileArg.getValue();
    } catch (TCLAP::ArgException& e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    try {
        const JXLImage image_in(filename_in);

        const SGEGBox& box = image_in.getBox();

        const bool has_multiple_spectral_sgeg = box.spectral_groups.size() > 1;

        for (const SGEGSpectralGroup& sg: box.spectral_groups) {
            std::string unique_filename_output = filename_out;

            if (has_multiple_spectral_sgeg) {
                // Edit the output path such as the output filename is unique
                std::string base, ext;
                Util::split_extension(filename_out.c_str(), base, ext);

                std::stringstream ss;
                ss << base << "_" << sg.root_name.data() << ext;

                unique_filename_output = ss.str();
            }

            if (sg.mins.size() != sg.maxs.size()) {
                throw std::runtime_error("Mins and maxs sizes don't match");
            }

            std::ofstream f_out(unique_filename_output);

            for (size_t i = 0; i < sg.mins.size(); i++) {
                f_out << sg.mins[i] << " " << sg.maxs[i] << std::endl;
            }
        }

    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

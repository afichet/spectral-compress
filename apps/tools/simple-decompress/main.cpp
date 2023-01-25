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
#include <chrono>
#include <fstream>

#include <EXRImage.h>
#include <JXLImage.h>

#include <tclap/CmdLine.h>


int main(int argc, char* argv[])
{
    std::string filename_in;
    std::string filename_out;

    bool save_timing = false;
    std::string filename_timing;

    try {
        TCLAP::CmdLine cmd("Utility to decompress JPEG-XL compressed with the simple-compress utility to OpenEXR.");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Specifies the JPEG-XL first input image.", true, "input.jxl", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the OpenEXR output image.", true, "output.exr", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<std::string> timingLogArg("l", "log", "Set file to save execution timing", false, "timing.txt", "path");

        cmd.add(timingLogArg);

        cmd.parse(argc, argv);

        filename_in  = inputFileArg.getValue();
        filename_out = outputFileArg.getValue();

        save_timing = timingLogArg.isSet();

        if (save_timing) {
            filename_timing = timingLogArg.getValue();
        }
    } catch (TCLAP::ArgException& e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    try {
        auto clock_start = std::chrono::steady_clock::now();

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

        auto clock_end = std::chrono::steady_clock::now();

        if (save_timing) {
            auto diff = clock_end - clock_start;

            std::ofstream logfile(filename_timing);
            logfile << "Total duration: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;
        }
    } catch (std::exception &e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

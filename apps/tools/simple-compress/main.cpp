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
#include <chrono>
#include <fstream>
#include <cstring>

#include <EXRImage.h>
#include <JXLImage.h>

#include <tclap/CmdLine.h>


int main(int argc, char* argv[])
{
    std::string filename_in;
    std::string filename_out;

    float framedistance = 0.f;

    bool save_timing = false;
    std::string filename_timing;

    try {
        TCLAP::CmdLine cmd("Utility to compress all OpenEXR framebuffer in JPEG-XL.");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Specifies the OpenEXR input image.", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the JPEG-XL output image.", true, "output.jxl", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<float> framedistanceArg("d", "framedistance", "Specifies the frame distance to use (compression parameter).", false, 0.f, "float");

        cmd.add(framedistanceArg);

        TCLAP::ValueArg<std::string> timingLogArg("l", "log", "Set file to save execution timing", false, "timing.txt", "path");

        cmd.add(timingLogArg);

        cmd.parse(argc, argv);

        filename_in  = inputFileArg.getValue();
        filename_out = outputFileArg.getValue();

        framedistance = framedistanceArg.getValue();

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

        const EXRImage exr_in(filename_in);
        JXLImage jxl_out(exr_in.width(), exr_in.height());
        SGEGBox box;

        box.exr_attributes = exr_in.getAttributesData();

        const std::vector<EXRFramebuffer*>& framebuffers = exr_in.getFramebuffersConst();

        for (const EXRFramebuffer* fb: framebuffers) {
            SGEGGrayGroup gg;
            const char* layer_name = fb->getName();

            std::pair<int, int> enc_bits;
            switch (fb->getPixelType()) {
                case Imf::PixelType::HALF:
                    enc_bits = std::make_pair(16, 5);
                    break;
                case Imf::PixelType::FLOAT:
                    enc_bits = std::make_pair(32, 8);
                    break;
                case Imf::PixelType::UINT:
                    // TODO: A rescaling is most probably necessary there
                    std::cout << "Warning: it is largely untested" << std::endl;
                    enc_bits = std::make_pair(32, 0);
                    break;
                default:
                    std::cerr << "Unknown OpenEXR type." << std::endl;
                    return 1;
            }

            gg.layer_index = jxl_out.appendFramebuffer(
                fb->getPixelDataConst(), 1,
                std::make_pair(32, 8), 1,
                framedistance,
                fb->getName());

            gg.layer_name.resize(std::strlen(layer_name) + 1);
            std::memcpy(gg.layer_name.data(), layer_name, gg.layer_name.size() * sizeof(char));

            box.gray_groups.push_back(gg);
        }

        jxl_out.setBox(box);

        jxl_out.write(filename_out);

        auto clock_end = std::chrono::steady_clock::now();

        if (save_timing) {
            auto diff = clock_end - clock_start;

            std::ofstream logfile(filename_timing);
            logfile << "Total duration: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;
        }
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

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
#include <fstream>

#include <tclap/CmdLine.h>

#include <EXRSpectralImage.h>
#include <Util.h>

int main(int argc, char* argv[])
{
    std::string filename_in, filename_out;
    size_t x, y;

    try {
        TCLAP::CmdLine cmd("Extract a spectrum from a given pixel location to a csv file.");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Specifies the spectral OpenEXR input image.", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies output CSV file.", true, "spectrum.csv", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<size_t> xCoordArg("x", "x_pos", "Sets the x coordinate to extract the spectrum from", true, 0, "coordinate");
        TCLAP::ValueArg<size_t> yCoordArg("y", "y_pos", "Sets the y coordinate to extract the spectrum from", true, 0, "coordinate");

        cmd.add(xCoordArg);
        cmd.add(yCoordArg);

        cmd.parse(argc, argv);

        filename_in  = inputFileArg.getValue();
        filename_out = outputFileArg.getValue();

        x = xCoordArg.getValue();
        y = yCoordArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    const EXRSpectralImage image_in(filename_in);

    const size_t width  = image_in.width();
    const size_t height = image_in.height();

    // Ensure the x and y coordinates are in range;
    if (x >= width || y >= height) {
        std::cerr << "Provided coordinates are out of range " << width << "x" << height << "px" << std::endl;
        return 1;
    }

    const std::vector<SpectralFramebuffer*>& spectral_buffers = image_in.getSpectralFramebuffersConst();
    const bool unique_output = spectral_buffers.size() == 1;

    for(const SpectralFramebuffer* fb: spectral_buffers) {
        std::string curr_filename_out;

        if (!unique_output) {
            std::string base, ext;
            Util::split_extension(filename_out.c_str(), base, ext);

            curr_filename_out = base + "_" + fb->root_name + ext;
        } else {
            curr_filename_out = filename_out;
        }

        const size_t n_bands = fb->wavelengths_nm.size();

        std::ofstream csv_file(curr_filename_out);

        for (size_t wl = 0; wl < n_bands; wl++) {
            csv_file << fb->wavelengths_nm[wl] << ", "
                     << fb->image_data[(y * width + x) * n_bands + wl]
                     << std::endl;
        }
    }

    return 0;
}

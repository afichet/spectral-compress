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

#include <string>
#include <cassert>
#include <cmath>

#include <EXRSpectralImage.h>
#include <Util.h>
#include <SpectrumConverter.h>

#include <tclap/CmdLine.h>
#include <lodepng.h>


inline float to_sRGB(float c)
{
    if (c <= 0.0031308f) {
        return 12.92f * c;
    } else {
        return 1.055 * std::pow(c, 1.f/2.4f) - 0.055;
    }
}


int main(int argc, char* argv[])
{
    std::string filename_input, filename_output;

    float exposure;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compares two Spectral OpenEXR images");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Specifies the EXR input image.", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the PNG output image.", true, "output.png", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<float> exposureArg("e", "exposure", "Sets the exposure value for the LDR conversion.", false, 0.f, "exposure");

        cmd.add(exposureArg);

        cmd.parse(argc, argv);

        filename_input  = inputFileArg.getValue();
        filename_output = outputFileArg.getValue();

        exposure  = exposureArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    EXRSpectralImage exr_in(filename_input);

    const size_t width  = exr_in.width();
    const size_t height = exr_in.height();

    const std::vector<SpectralFramebuffer*>& spectral_framebuffers = exr_in.getSpectralFramebuffers();

    const bool has_multiple_framebuffers = spectral_framebuffers.size() > 1;

    SpectrumConverter reflective_converter(false);
    SpectrumConverter emissive_converter(true);

    for (const SpectralFramebuffer* fb: spectral_framebuffers) {
        std::string unique_filename_output = filename_output;

        if (has_multiple_framebuffers) {
            // Edit the output path such as the output filename is unique
            std::string base, ext;
            Util::split_extension(filename_output.c_str(), base, ext);

            std::stringstream ss;
            ss << base << "_" << fb->root_name << ext;

            unique_filename_output = ss.str();
        }

        // Do the conversion
        std::vector<float> rgb_image;
        bool hasRGB = false;

        if (Util::ends_with(fb->root_name, "S0")) {
            emissive_converter.spectralImageToRGB(
                fb->wavelengths_nm,
                fb->image_data,
                width, height,
                rgb_image
            );

            hasRGB = true;
        } else if (Util::ends_with(fb->root_name, "T")) {
            reflective_converter.spectralImageToRGB(
                fb->wavelengths_nm,
                fb->image_data,
                width, height,
                rgb_image
            );

            hasRGB = true;
        }

        if (hasRGB) {
            assert(rgb_image.size() == 3 * width * height);

            const float exposure_v = std::exp2(exposure);

            // Apply the exposure compensation & sRGB gamma
            #pragma omp parallel for
            for (size_t i = 0; i < 3 * width * height; i++) {
                rgb_image[i] *= exposure_v;
                rgb_image[i] = to_sRGB(rgb_image[i]);
            }

            // Creates an RGB visualisation
            std::vector<uint8_t> framebuffer_rgba(4 * width * height);

            #pragma omp parallel for
            for (size_t px = 0; px < width * height; px++) {
                for (size_t c = 0; c < 3; c++) {
                    framebuffer_rgba[4 * px + c] =
                        std::round(255.f * std::max(0.f, std::min(1.f, rgb_image[3 * px + c])));
                }

                framebuffer_rgba[4 * px + 3] = 255;
            }

            // Saves the file
            lodepng::encode(unique_filename_output, framebuffer_rgba, width, height);
        }
    }

    return 0;
}

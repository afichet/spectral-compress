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
#include <cassert>
#include <algorithm>
#include <chrono>
#include <fstream>

#include <JXLImage.h>
#include <EXRSpectralImage.h>
#include <Util.h>
#include <moments_image.h>

#include <tclap/CmdLine.h>


void decompress_spectral_framebuffer(
    const SGEGSpectralGroup& sg,
    const std::vector<std::vector<float>>& compressed_moments,
    std::vector<float>& spectral_framebuffer)
{
    const size_t n_pixels = compressed_moments[0].size();
    const size_t n_bands  = sg.wavelengths.size();

    size_t n_moments = 0;

    assert(sg.layer_indices.size() == (sg.mins.size() + 1));
    assert(sg.layer_indices.size() == (sg.maxs.size() + 1));
    n_moments = compressed_moments.size();

    std::vector<double> wavelengths_d;
    std::vector<double> mins_d, maxs_d;
    std::vector<double> spectral_framebuffer_d;
    std::vector<double> compressed_moments_d(n_moments * n_pixels);

    Util::cast_vector(sg.wavelengths, wavelengths_d);
    Util::cast_vector(sg.mins, mins_d);
    Util::cast_vector(sg.maxs, maxs_d);

    for (size_t px = 0; px < n_pixels; px++) {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments_d[px * n_moments + m] = compressed_moments[m][px];
        }
    }

    double timing;

    decompress_spectral_image(
        sg.method,
        wavelengths_d,
        compressed_moments_d,
        mins_d, maxs_d,
        n_pixels,
        n_moments,
        spectral_framebuffer_d,
        timing
    );

    const std::string root_name(sg.root_name.data());

    if (Util::ends_with(root_name, "S0") || Util::ends_with(root_name, "T")) {
        spectral_framebuffer.resize(n_pixels * n_bands);

        // Ensure positive values
        #pragma omp parallel for
        for (size_t i = 0; i < n_pixels * n_bands; i++) {
            spectral_framebuffer[i] = (float)std::max(spectral_framebuffer_d[i], 0.);
        }
    } else {
        Util::cast_vector(spectral_framebuffer_d, spectral_framebuffer);
    }
}


Imf::PixelType quantization_to_exr(size_t n_bits, size_t n_exponent_bits) {
    if (n_exponent_bits > 0) {
        if (n_bits > 16) {
            return Imf::PixelType::FLOAT;
        } else {
            return Imf::PixelType::HALF;
        }
    } else {
        return Imf::PixelType::UINT;
    }
}


int main(int argc, char* argv[])
{
    std::string filename_in;
    std::string filename_out;

    bool save_timing = false;
    std::string filename_timing;

    try {
        TCLAP::CmdLine cmd("Utility to decompress JPEG-XL compressed with the compress utility to OpenEXR.");

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

        const JXLImage jxl_image(filename_in);

        const SGEGBox box   = jxl_image.getBox();
        const size_t width  = jxl_image.width();
        const size_t height = jxl_image.height();

        EXRSpectralImage exr_out(width, height);

        exr_out.setAttributesData(box.exr_attributes);

        for (const SGEGSpectralGroup& sg: box.spectral_groups) {
            std::string root_name = sg.root_name.data();
            const size_t n_moments = sg.layer_indices.size();

            assert(sg.mins.size() == sg.maxs.size());

            std::vector<std::vector<float>> moments(n_moments);
            std::vector<float> spectral_framebuffer;

            for (size_t m = 0; m < n_moments; m++) {
                moments[m] = jxl_image.getFramebufferDataConst(sg.layer_indices[m]);
            }

            decompress_spectral_framebuffer(
                sg,
                moments,
                spectral_framebuffer
            );

            const JXLFramebuffer* main_fb = jxl_image.getFramebuffer(sg.layer_indices[0]);
            const size_t n_bits           = main_fb->getBitsPerSample();
            const size_t n_exponent_bits  = main_fb->getExponentBitsPerSample();
            const Imf::PixelType pixel_type = quantization_to_exr(n_bits, n_exponent_bits);

            exr_out.appendSpectralFramebuffer(
                sg.wavelengths,
                spectral_framebuffer,
                root_name,
                pixel_type
            );
        }

        for (const SGEGGrayGroup& gg: box.gray_groups) {
            const JXLFramebuffer* main_fb = jxl_image.getFramebuffer(gg.layer_index);
            const size_t n_bits           = main_fb->getBitsPerSample();
            const size_t n_exponent_bits  = main_fb->getExponentBitsPerSample();
            const Imf::PixelType pixel_type = quantization_to_exr(n_bits, n_exponent_bits);

            exr_out.appendExtraFramebuffer(
                jxl_image.getFramebufferDataConst(gg.layer_index),
                gg.layer_name.data(),
                pixel_type
            );
        }

        exr_out.write(filename_out);

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

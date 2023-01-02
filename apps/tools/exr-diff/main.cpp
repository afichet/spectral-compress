/**
 * Copyright 2022 Alban Fichet
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
#include <sstream>
#include <cstdint>
#include <cstddef>
#include <cassert>

#include <EXRSpectralImage.h>
#include <Util.h>

#include <tclap/CmdLine.h>

#include "colormapdata.h"
#include <lodepng.h>


bool compare_spectral_images(
    const SpectralFramebuffer* img_a,
    const SpectralFramebuffer* img_b,
    uint32_t width, uint32_t height,
    std::vector<float>& diff_image,
    float& min_val, float& max_val,
    float& avg_err)
{
    // Check if dimensions match
    if (img_a->wavelengths_nm.size() != img_b->wavelengths_nm.size()) {
        return false;
    }

    const size_t n_bands = img_a->wavelengths_nm.size();

    // Compute image diff
    const std::vector<float>& data_a = img_a->image_data;
    const std::vector<float>& data_b = img_b->image_data;

    assert(data_a.size() == data_b.size());
    diff_image.resize(width * height);

    #pragma omp parallel for
    for (size_t px = 0; px < width * height; px++) {
        diff_image[px] = 0;

        float diff = 0;

        for (size_t b = 0; b < n_bands; b++) {
            float d = data_a[px * n_bands + b] - data_b[px * n_bands + b];

            diff += d*d;
        }

        diff_image[px] = std::sqrt(diff) / (float)n_bands;
    }

    // Compute min max (cannot be multi-threaded)
    min_val = diff_image[0];
    max_val = diff_image[0];

    for (size_t i = 1; i < width * height; i++) {
        min_val = std::min(min_val, diff_image[i]);
        max_val = std::max(max_val, diff_image[i]);
    }

    avg_err = 0;

    for (size_t px = 0; px < width * height; px++) {
        avg_err += diff_image[px];
    }

    avg_err /= (float)(width * height);

    return true;
}


void diff_to_rgba(
    const std::vector<float>& diff_image,
    float lower_bound, float higher_bound,
    const float lut[],
    size_t lut_size,
    std::vector<uint8_t>& rgba_image)
{
    assert(lower_bound < higher_bound);

    rgba_image.resize(4 * diff_image.size());

    for (size_t i = 0; i < diff_image.size(); i++) {
        const float v = std::max(0.f, std::min(1.f, (diff_image[i] - lower_bound) / (higher_bound - lower_bound)));
        const size_t idx_lut = (size_t)std::round(v * lut_size);

        rgba_image[4 * i + 0] = std::round(255.f * v * lut[3 * idx_lut + 0]);
        rgba_image[4 * i + 1] = std::round(255.f * v * lut[3 * idx_lut + 1]);
        rgba_image[4 * i + 2] = std::round(255.f * v * lut[3 * idx_lut + 2]);
        rgba_image[4 * i + 3] = 255;
    }
}


int main(int argc, char* argv[])
{
    std::string filename_a, filename_b;
    std::string filename_output;

    float custom_lower_bound, custom_upper_bound;

    bool custom_lower_bound_is_set = false;
    bool custom_upper_bound_is_set = false;

    std::string error_output;

    bool error_output_is_set = false;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compares two Spectral OpenEXR images");

        TCLAP::UnlabeledValueArg<std::string> inputAFileArg("FileA", "Specifies reference image (input).", true, "reference.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> inputBFileArg("FileB", "Specifies image to compare with (input).", true, "comparison.jxl", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the image to write into (output).", true, "output.png", "path");

        cmd.add(inputAFileArg);
        cmd.add(inputBFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<float> lowerBoundArg("l", "lower", "Sets the lower bound for the colormap.", false, 0.f, "min");
        TCLAP::ValueArg<float> upperBoundArg("u", "upper", "Sets the upper bouind for the colormap.", false, 0.f, "max");

        cmd.add(lowerBoundArg);
        cmd.add(upperBoundArg);

        TCLAP::ValueArg<std::string> errorFileArg("e", "error", "Sets a file where to write error log.", false, "error.bin", "path");

        cmd.add(errorFileArg);

        cmd.parse(argc, argv);

        filename_a = inputAFileArg.getValue();
        filename_b = inputBFileArg.getValue();
        filename_output = outputFileArg.getValue();

        custom_lower_bound_is_set = lowerBoundArg.isSet();
        custom_upper_bound_is_set = upperBoundArg.isSet();

        custom_lower_bound  = lowerBoundArg.getValue();
        custom_upper_bound = upperBoundArg.getValue();

        error_output_is_set = errorFileArg.isSet();
        error_output = errorFileArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    FILE* f_err = NULL;

    if (error_output_is_set) {
        f_err = fopen(error_output.c_str(), "wb");
    }

    EXRSpectralImage exr_in_a(filename_a);
    EXRSpectralImage exr_in_b(filename_b);

    const size_t width  = exr_in_a.width();
    const size_t height = exr_in_a.height();

    if (exr_in_b.width() != width || exr_in_b.height() != height) {
        std::cerr << "Files dimensions do not match!" << std::endl;
        return 1;
    }

    const std::vector<SpectralFramebuffer*>& spectral_framebuffers_a = exr_in_a.getSpectralFramebuffers();
    const std::vector<SpectralFramebuffer*>& spectral_framebuffers_b = exr_in_b.getSpectralFramebuffers();

    const bool has_multiple_framebuffers = spectral_framebuffers_a.size() > 1;

    for (const SpectralFramebuffer* fb_a: spectral_framebuffers_a) {
        // Now, look for the same rootname in the other file
        for (const SpectralFramebuffer* fb_b: spectral_framebuffers_b) {
            if (fb_a->root_name == fb_b->root_name) {
                std::string unique_filename_output = filename_output;

                if (has_multiple_framebuffers) {
                    // Edit the output path such as the output filename is unique
                    std::string base, ext;
                    Util::split_extension(filename_output.c_str(), base, ext);

                    std::stringstream ss;
                    ss << base << "_" << fb_a->root_name << ext;

                    unique_filename_output = ss.str();
                }

                // Do the comparison
                std::vector<float> framebuffer_error;
                float min_err, max_err, avg_err;

                compare_spectral_images(
                    fb_a, fb_b,
                    width, height,
                    framebuffer_error,
                    min_err, max_err,
                    avg_err
                );

                if (error_output_is_set) {
                    fwrite(&avg_err, sizeof(float), 1, f_err);
                }
                // std::cout << "Error: [" << min_err << ", " << max_err << "]" << std::endl;

                float lower, upper;

                if (custom_upper_bound_is_set) {
                    lower = custom_lower_bound;
                } else {
                    lower = 0;
                }

                if (custom_lower_bound_is_set) {
                    upper = custom_upper_bound;
                } else {
                    upper = max_err;
                }

                // Creates an RGB visualisation
                std::vector<uint8_t> framebuffer_rgba;

                diff_to_rgba(
                    framebuffer_error,
                    lower, upper,
                    turbo_colormap_data,
                    sizeof(turbo_colormap_data) / (3 * sizeof(float)),
                    framebuffer_rgba
                );

                // Saves the file
                lodepng::encode(unique_filename_output, framebuffer_rgba, width, height);

                break;
            }
        }
    }

    if (error_output_is_set) {
        fclose(f_err);
    }

    return 0;
}

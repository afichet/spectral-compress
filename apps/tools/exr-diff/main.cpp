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
#include <sstream>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <stdexcept>
#include <algorithm>

#include <EXRSpectralImage.h>
#include <Util.h>

#include <tclap/CmdLine.h>

#include "colormapdata.h"
#include <lodepng.h>


void rmse_spectral_images(
    const SpectralFramebuffer* img_a,
    const SpectralFramebuffer* img_b,
    size_t width, size_t height,
    std::vector<double>& rmse_image,
    double& min_val, double& max_val,
    double& avg_err,
    double& five_higher_percentile)
{
    // Check if dimensions match
    if (img_a->wavelengths_nm.size() != img_b->wavelengths_nm.size()) {
        throw std::runtime_error("dimensions mismatch");
    }

    const size_t n_bands = img_a->wavelengths_nm.size();

    // Compute image diff
    const std::vector<float>& data_a = img_a->image_data;
    const std::vector<float>& data_b = img_b->image_data;

    assert(data_a.size() == data_b.size());
    rmse_image.resize(width * height);

    #pragma omp parallel for
    for (size_t px = 0; px < width * height; px++) {
        rmse_image[px] = 0;

        double diff = 0;

        for (size_t b = 0; b < n_bands; b++) {
            double d =
                (double)data_a[px * n_bands + b]
              - (double)data_b[px * n_bands + b];

            diff += d*d;
        }

        rmse_image[px] = std::sqrt(diff / (double)n_bands);
    }

    // Compute min max (cannot be multi-threaded)
    min_val = rmse_image[0];
    max_val = rmse_image[0];

    for (size_t i = 1; i < width * height; i++) {
        min_val = std::min(min_val, rmse_image[i]);
        max_val = std::max(max_val, rmse_image[i]);
    }

    avg_err = 0;

    for (size_t px = 0; px < width * height; px++) {
        avg_err += rmse_image[px];
    }

    avg_err /= (double)(width * height);

    // Calculate 5% higher
    std::vector<double> rmses(rmse_image);
    std::sort(rmses.begin(), rmses.end());

    size_t idx = std::round((float)rmses.size() * 95.f / 100.f);
    five_higher_percentile = rmses[idx];

    std::cout << "average rmse / pixel: " << avg_err << std::endl;
    std::cout << "    max rmse / pixel: " << max_val << std::endl;
    // std::cout << "         rmse global: " << Util::rmse_images(data_a, data_b, width * height, n_bands) << std::endl;
    // std::cout << "        rrmse global: " << Util::rrmse_images(data_a, data_b, width * height, n_bands) << std::endl;
    // std::cout << "          max / band: " << Util::max_error_images(data_a, data_b, width * height, n_bands) << std::endl;
    std::cout << "    5 percent higher: " << five_higher_percentile << std::endl;
}


template<typename T>
void diff_to_rgba(
    const std::vector<T>& diff_image,
    float lower_bound, float upper_bound,
    const float lut[],
    size_t lut_size,
    std::vector<uint8_t>& rgba_image)
{
    assert(lower_bound < upper_bound);

    rgba_image.resize(4 * diff_image.size());

    for (size_t i = 0; i < diff_image.size(); i++) {
        const T v = Util::clamp((diff_image[i] - lower_bound) / (upper_bound - lower_bound), (T)0, (T)1);
        const size_t idx_lut = (size_t)std::round(v * (lut_size - 1));

        assert(idx_lut < lut_size);

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

    double custom_lower_bound, custom_upper_bound;

    bool custom_lower_bound_is_set = false;
    bool custom_upper_bound_is_set = false;

    std::string error_output;

    bool error_output_is_set = false;

    bool generate_scale = false;
    uint32_t scale_width, scale_height;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compares two Spectral OpenEXR images showing an heatmap of the RMSE per pixel");

        TCLAP::UnlabeledValueArg<std::string> inputAFileArg("FileA", "Specifies reference image (input).", true, "reference.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> inputBFileArg("FileB", "Specifies image to compare with (input).", true, "comparison.jxl", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the image to write into (output).", true, "output.png", "path");

        cmd.add(inputAFileArg);
        cmd.add(inputBFileArg);
        cmd.add(outputFileArg);

        TCLAP::ValueArg<double> lowerBoundArg("l", "lower", "Sets the lower bound for the colormap.", false, 0.f, "min");
        TCLAP::ValueArg<double> upperBoundArg("u", "upper", "Sets the upper bouind for the colormap.", false, 1.f, "max");

        cmd.add(lowerBoundArg);
        cmd.add(upperBoundArg);

        TCLAP::ValueArg<std::string> errorFileArg("e", "error", "Sets a file where to write error log.", false, "error.bin", "path");

        cmd.add(errorFileArg);

        TCLAP::SwitchArg generateScaleArg("g", "generate", "Generate a scale instead of comparing the images.", false);
        TCLAP::ValueArg<uint32_t> generateScaleWArg("x", "width", "Scale width. Works with -g parameter", false, 30, "width");
        TCLAP::ValueArg<uint32_t> generateScaleHArg("y", "height", "Scale height. Works with -g parameter", false, 300, "height");
        // FIXME: this is not the nicest way of doing things, constraints shall
        // be added to have mutually exclusive args, namely:
        // inputAFileArg, inputBFileArg, lowerBoundArg, upperBoundArg, errorFileArg
        // vs. generateScaleArg
        cmd.add(generateScaleArg);
        cmd.add(generateScaleWArg);
        cmd.add(generateScaleHArg);

        cmd.parse(argc, argv);

        filename_a = inputAFileArg.getValue();
        filename_b = inputBFileArg.getValue();
        filename_output = outputFileArg.getValue();

        custom_lower_bound_is_set = lowerBoundArg.isSet();
        custom_upper_bound_is_set = upperBoundArg.isSet();

        custom_lower_bound = lowerBoundArg.getValue();
        custom_upper_bound = upperBoundArg.getValue();

        error_output_is_set = errorFileArg.isSet();
        error_output = errorFileArg.getValue();

        generate_scale = generateScaleArg.getValue();

        scale_width = generateScaleWArg.getValue();
        scale_height = generateScaleHArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arguemnt " << e.argId() << std::endl;
        return 1;
    }

    if (generate_scale) {
        // Just generates a scale
        std::vector<float> scale_values(scale_width * scale_height);
        std::vector<uint8_t> scale_rgba(4 * scale_width * scale_height);

        for (uint32_t y = 0; y < scale_height; y++) {
            const float v = float(scale_height - y - 1) / float(scale_height - 1);

            for (uint32_t x = 0; x < scale_width; x++) {
                scale_values[y * scale_width + x] = v;
            }
        }

        diff_to_rgba(
            scale_values,
            0.f, 1.f,
            turbo_colormap_data,
            sizeof(turbo_colormap_data) / (3 * sizeof(float)),
            scale_rgba
        );

        lodepng::encode(filename_output, scale_rgba, scale_width, scale_height);
    } else {
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
                    std::vector<double> framebuffer_error;
                    double min_err, max_err, avg_err, five_higher_percentile;

                    std::cout << "Statistics for " << fb_a->root_name << std::endl;

                    rmse_spectral_images(
                        fb_a, fb_b,
                        width, height,
                        framebuffer_error,
                        min_err, max_err,
                        avg_err,
                        five_higher_percentile
                    );

                    if (error_output_is_set) {
                        fwrite(&avg_err, sizeof(double), 1, f_err);
                        fwrite(&five_higher_percentile, sizeof(double), 1, f_err);
                    }

                    float lower, upper;

                    if (custom_lower_bound_is_set) {
                        lower = custom_lower_bound;
                    } else {
                        lower = 0;
                    }

                    if (custom_upper_bound_is_set) {
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
    }

    return 0;
}

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

#include "stats.h"
#include "cmd_constraints.h"

#include <JXLImage.h>
#include <EXRSpectralImage.h>

#include <moments.h>
#include <moments_image.h>
#include <curve_quantization.h>
#include <curve_compression.h>

#include <tclap/CmdLine.h>

#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <chrono>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstring>

/**
 * TODO:
 * - Add an argument
 *    - to control spatial downsampling (needs as well serious work in the JXL class)
 * - RGB layers are only computed on the root (maybe a desirable behaviour though...)
 */


void compress_spectral_framebuffer(
    SpectralCompressionType method,
    // Input data
    const SpectralFramebuffer* framebuffer,
    uint32_t width, uint32_t height,
    // Quantization
    std::pair<int, int> n_bits_dc,
    std::pair<int, int> n_bits_ac1,
    bool uses_constant_quantization,
    std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    // Compression
    float compression_dc,
    float compression_ac1,
    bool uses_constant_compression,
    std::vector<float>& compression_curve,
    int effort,
    // Downsampling
    uint32_t downsampling_factor_ac1,
    std::vector<uint32_t>& downsampling_factor_curve,
    // Result of compression
    std::vector<std::vector<float>>& compressed_moments,
    SGEGSpectralGroup& sg,
    // Logging
    bool log,
    std::stringstream& log_stream,
    stats_data &quantization_error,
    stats_data &compression_error)
{
    const uint32_t n_moments = framebuffer->wavelengths_nm.size();
    const uint32_t n_pixels = width * height;

    assert(framebuffer->image_data.size() == n_moments * n_pixels);

    std::vector<double> compressed_moments_d;
    std::vector<uint8_t> relative_scales;
    std::vector<double> mins_d, maxs_d;
    double global_min_d, global_max_d;

    // TODO: this shall be revisited
    // TODO: this is invalid for polarization layers that can be negative
    // Condition buffers to ensure everything is in order
    std::vector<double> wavelengths(framebuffer->wavelengths_nm.size());

    #pragma omp parallel for
    for (size_t i = 0; i < wavelengths.size(); i++) {
        wavelengths[i] = framebuffer->wavelengths_nm[i];
    }

    std::vector<double> spectral_image(framebuffer->image_data.size());

    const char back_char = framebuffer->root_name.back();

    // Check if we're dealing with a polarisation component
    if (back_char == '1' || back_char == '2' || back_char == '3') {
        // TODO
        std::cout << "Warning: we are dealing with a polarisation layer, this is not yet supported and the result will be wrong." << std::endl;
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < spectral_image.size(); i++) {
            double v = framebuffer->image_data[i];

            if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
                v = 1e-8;
            }

            spectral_image[i] = v;
        }
    }

    double quantization_curve_timing(0),
           compression_curve_timing(0),
           compression_timing(0);

    compute_quantization_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        n_bits_dc,
        n_bits_ac1,
        uses_constant_quantization,
        quantization_curve,
        normalize_moments,
        quantization_curve_timing
    );

    // TODO: We don't automatically compute the downsampling curves yet
    downsampling_factor_curve.resize(n_moments);
    downsampling_factor_curve[0] = 1;

    for (size_t i = 1; i < n_moments; i++) {
        downsampling_factor_curve[i] = downsampling_factor_ac1;
    }

    compute_compression_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        normalize_moments,
        downsampling_factor_curve,
        compression_dc, compression_ac1,
        uses_constant_compression,
        compression_curve,
        effort,
        compression_curve_timing
    );

    compress_spectral_image(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        compressed_moments_d,
        mins_d, maxs_d,
        normalize_moments,
        relative_scales,
        global_min_d, global_max_d,
        compression_timing
    );

    // Compute errors (may be done twice if the quantization / compression
    // profile is built at runtime)

    quantization_error = stats_for_quantization_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        compressed_moments_d, mins_d, maxs_d,
        relative_scales, global_min_d, global_max_d
    );

    compression_error = stats_for_compression_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        downsampling_factor_curve,
        compression_curve,
        compressed_moments_d, mins_d, maxs_d,
        relative_scales, global_min_d, global_max_d,
        effort
    );

    // Give back data to the caller

    compressed_moments.resize(n_moments);

    for (size_t m = 0; m < n_moments; m++) {
        compressed_moments[m].resize(n_pixels);

        for (size_t px = 0; px < n_pixels; px++) {
            compressed_moments[m][px] = (float)compressed_moments_d[n_moments * px + m];
        }
    }

    sg.method = method;

    sg.mins.resize(n_moments - 1);
    sg.maxs.resize(n_moments - 1);

    for (size_t m = 0; m < n_moments - 1; m++) {
        sg.mins[m] = (float)mins_d[m];
        sg.maxs[m] = (float)maxs_d[m];
    }

    // Copy the local scaling
    if (method == UPPERBOUND || method == TWOBOUNDS) {
        assert(relative_scales.size() == n_pixels);

        // Extend by 1 for the extra scaling framebuffer
        quantization_curve.push_back(std::make_pair(8, 0));
        compression_curve.push_back(compression_dc);
        downsampling_factor_curve.push_back(1);

        compressed_moments.resize(n_moments + 1);
        compressed_moments[n_moments].resize(n_pixels);

        for (size_t px = 0; px < n_pixels; px++) {
            compressed_moments[n_moments][px] = (float)((double)relative_scales[px] / (double)std::numeric_limits<uint8_t>::max());
        }

        sg.global_min = global_min_d;
        sg.global_max = global_max_d;
    }

    if (log) {
        switch(sg.method) {
            case LINEAR:
                log_stream << "Method: linear" << std::endl;
                break;
            case LINAVG:
                log_stream << "Method: linavg" << std::endl;
                break;
            case BOUNDED:
                log_stream << "Method: bounded" << std::endl;
                break;
            case UNBOUNDED:
                log_stream << "Method: unbounded" << std::endl;
                break;
            case UNBOUNDED_TO_BOUNDED:
                log_stream << "Method: unbounded_to_bounded" << std::endl;
                break;
            case UPPERBOUND:
                log_stream << "Method: upperbound" << std::endl;
                break;
            case TWOBOUNDS:
                log_stream << "Method: twobounds" << std::endl;
                break;
        }
        log_stream << "Quantization curve:" << std::endl;
        for (const std::pair<int, int>& q : quantization_curve) {
            log_stream << q.first << " ";
        }
        log_stream << std::endl;
        for (const std::pair<int, int>& q : quantization_curve) {
            log_stream << q.second << " ";
        }
        log_stream << std::endl;

        log_stream << "rmse: "  << quantization_error.rmse_error << std::endl;
        log_stream << "rrmse: " << quantization_error.rrmse_error << std::endl;
        log_stream << "avg: "   << quantization_error.avg_error << std::endl;
        log_stream << "max: "   << quantization_error.max_error << std::endl;

        log_stream << "Compression curve:" << std::endl;
        for (const float& c : compression_curve) {
            log_stream << c << " ";
        }
        log_stream << std::endl;
        log_stream << "rmse: "  << compression_error.rmse_error << std::endl;
        log_stream << "rrmse: " << compression_error.rrmse_error << std::endl;
        log_stream << "avg: "   << compression_error.avg_error << std::endl;
        log_stream << "max: "   << compression_error.max_error << std::endl;

        log_stream << "Timings:" << std::endl;
        log_stream << "Quantization curve: " << quantization_curve_timing << " ms" << std::endl;
        log_stream << "Compression curve:  " << compression_curve_timing << " ms" << std::endl;
        log_stream << "Compression:        " << compression_timing << " ms" << std::endl;
    }
}


void quantization_from_exr(Imf::PixelType type, std::pair<int, int>& n_bits) {
    switch (type) {
        case Imf::PixelType::UINT:
            n_bits.first = 32;
            n_bits.second = 0;
            break;

        case Imf::PixelType::HALF:
            n_bits.first = 16;
            n_bits.second = 5;
            break;

        case Imf::PixelType::FLOAT:
            n_bits.first = 32;
            n_bits.second = 8;
            break;

        default:
            throw std::runtime_error("Unknown pixel type");
            break;
    }
}


int main(int argc, char *argv[])
{
    std::string filename_in, filename_out, filename_dump;
    bool log_is_active = false;

    bool has_txt_log;
    bool has_bin_log;
    bool has_dump_file;

    std::string log_filepath, binarylog_filepath;
    std::stringstream log_content;

    float compression_dc = .1f;
    float compression_ac1 = .1f;
    bool use_flat_compression = false;
    int compression_effort = 7;

    std::pair<int, int> n_bits_ac1;
    bool use_flat_quantization = false;
    bool normalize_moments = true;

    int downsampling_factor_ac1 = 1;

    SpectralCompressionType method = TWOBOUNDS;

    std::chrono::time_point<std::chrono::steady_clock>  clock_start, clock_end;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compress a spectral image");

        // Input / output
        TCLAP::UnlabeledValueArg<std::string> inputFileArg("Input", "Specifies the Spectral OpenEXR file to compress from (input).", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Specifies the JPEG XL file to compress spectral data to (output).", true, "output.jxl", "path");
        TCLAP::ValueArg<std::string> logFileArg("l", "log", "Specifies a file to log timings and quantization and compression curves into,", false, "log.txt", "path");
        TCLAP::ValueArg<std::string> binarylogFileArg("k", "binary_log", "Specifies a file to log timings into in binary form,", false, "log.bin", "path");
        TCLAP::ValueArg<std::string> dumpFileArg("d", "dump", "Saves a dump file before applying JXL compression", false, "dump.bin", "path");
        cmd.add(inputFileArg);
        cmd.add(outputFileArg);
        cmd.add(logFileArg);
        cmd.add(binarylogFileArg);
        cmd.add(dumpFileArg);

        // Compresion tweaking
        FrameDistanceConstraint frameDistanceConstraint;
        CompressionEffortConstraint compressionEffortConstraint;

        TCLAP::ValueArg<float> frameDistanceDCArg("a", "frame_distance_dc", "Sets the distance level for lossy compression (compression rate) on the DC component.", false, .1f, &frameDistanceConstraint);
        TCLAP::ValueArg<float> frameDistanceACArg("b", "frame_distance_ac", "Sets the distance level for lossy compression (compression rate) on the first AC component. The program uses the same compression ratio for the remaining components when `--c_flat` is set. Otherwise, the program generates a compression curve starting with the distance parameter for the first AC component using the provided value based on the image data (can be slow).", false, .1f, &frameDistanceConstraint);
        TCLAP::SwitchArg useFlatCompressionArg("c", "c_flat", "Sets a flat compression curve. All AC components use the same distance level as the first AC component set by `--frame_distance_ac` while the DC components uses the distance level provided by `--frame_distance_dc` parameter.");
        TCLAP::ValueArg<int> compressionEffortArg("e", "effort", "Sets the compression effort: 1 = fast, 9 = slow.", false, 7, &compressionEffortConstraint);
        cmd.add(frameDistanceDCArg);
        cmd.add(frameDistanceACArg);
        cmd.add(useFlatCompressionArg);
        cmd.add(compressionEffortArg);

        // Quantization tweaking
        TCLAP::ValueArg<int> quantizationMainStartArg("q", "quantization_bits", "Sets the starting number of bits for quantizing the first AC component. The program use the same number of bits for the remaining components when `--q_flat` is set. Otherwise, the program generates a custom quantization curve starting with the number of bits for the first AC component using the provided value based on the image data (can be slow).", false, 8, "integer");
        TCLAP::ValueArg<int> quantizationExpoStartArg("r", "quantization_exp", "Sets the starting number of bits for quantizing the first AC component. The program use the same number of bits for the remaining components when `--q_flat` is set. Otherwise, the program generates a custom quantization curve starting with the number of bits for the first AC component using the provided value based on the image data (can be slow).", false, 0, "integer");
        TCLAP::SwitchArg useFlatQuantizationArg("u", "q_flat", "Sets a flat quantization curve. The DC component uses the same quantization as the one used for all spectral data in the OpenEXR file while the remaining AC components use the number of bits provided in `--quantization parameter`.");
        TCLAP::SwitchArg normalizeMomentsArg("z", "kill_normalize", "Do not normalize the moment images. May cause problem when using integer framebuffer for storing moments.", false);
        cmd.add(quantizationMainStartArg);
        cmd.add(quantizationExpoStartArg);
        cmd.add(useFlatQuantizationArg);
        cmd.add(normalizeMomentsArg);

        // Downsampling tweaking
        DownsamplingFactorConstraint downsamplingFactorConstraint;

        TCLAP::ValueArg<int> downsamplingFactorArg("s", "downsampling", "Sets the downsampling ratio to apply to AC components.", false, 1, &downsamplingFactorConstraint);
        cmd.add(downsamplingFactorArg);

        // Moment storage method tweaking
        std::vector<std::string> allowedCompressionMethods;
        allowedCompressionMethods.push_back("linear");
        allowedCompressionMethods.push_back("linavg");
        allowedCompressionMethods.push_back("bounded");
        allowedCompressionMethods.push_back("unbounded");
        allowedCompressionMethods.push_back("unbounded_to_bounded");
        allowedCompressionMethods.push_back("upperbound");
        allowedCompressionMethods.push_back("twobounds");
        TCLAP::ValuesConstraint<std::string> allowedCompressionMethodsVals(allowedCompressionMethods);
        TCLAP::ValueArg<std::string> momentCompressionMethodArg("m", "method", "Representation of moments to use.", false, "twobounds", &allowedCompressionMethodsVals);
        cmd.add(momentCompressionMethodArg);

        cmd.parse(argc, argv);

        filename_in     = inputFileArg.getValue();
        filename_out    = outputFileArg.getValue();
        has_txt_log     = logFileArg.isSet();
        has_bin_log     = binarylogFileArg.isSet();
        log_is_active   = logFileArg.isSet() || binarylogFileArg.isSet();
        has_dump_file   = dumpFileArg.isSet();

        if (has_txt_log) {
            log_filepath = logFileArg.getValue();
        }

        if (has_bin_log) {
            binarylog_filepath = binarylogFileArg.getValue();
        }

        if (has_dump_file) {
            filename_dump = dumpFileArg.getValue();
        }

        compression_dc       = frameDistanceDCArg.getValue();
        compression_ac1      = frameDistanceACArg.getValue();
        use_flat_compression = useFlatCompressionArg.getValue();
        compression_effort   = compressionEffortArg.getValue();

        n_bits_ac1.first      = quantizationMainStartArg.getValue();
        n_bits_ac1.second     = quantizationExpoStartArg.getValue();
        use_flat_quantization = useFlatQuantizationArg.getValue();
        // normalize_moments     = !normalizeMomentsArg.getValue();

        downsampling_factor_ac1 = downsamplingFactorArg.getValue();

        if (momentCompressionMethodArg.getValue() == "linear") {
            method = LINEAR;
        } else if (momentCompressionMethodArg.getValue() == "linavg") {
            method = LINAVG;
        } else if (momentCompressionMethodArg.getValue() == "bounded") {
            method = BOUNDED;
        } else if (momentCompressionMethodArg.getValue() == "unbounded") {
            method = UNBOUNDED;
        } else if (momentCompressionMethodArg.getValue() == "unbounded_to_bounded") {
            method = UNBOUNDED_TO_BOUNDED;
        } else if (momentCompressionMethodArg.getValue() == "upperbound") {
            method = UPPERBOUND;
        } else if (momentCompressionMethodArg.getValue() == "twobounds") {
            method = TWOBOUNDS;
        }
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;

        return 1;
    }

    // Timing setup
    clock_start = std::chrono::steady_clock::now();

    // Load input file
    const EXRSpectralImage exr_in(filename_in);

    const std::vector<SpectralFramebuffer*>& spectral_framebuffers = exr_in.getSpectralFramebuffersConst();
    const std::vector<GreyFramebuffer*>& extra_framebuffers        = exr_in.getExtraFramebuffersConst();

    // Create output file
    JXLImage jxl_out(exr_in.width(), exr_in.height());
    SGEGBox box;

    // Errors
    std::vector<stats_data> quantization_errors, compression_errors;

    box.exr_attributes = exr_in.getAttributesData();

    // Run compression for each spectral group
    for (const SpectralFramebuffer* fb: spectral_framebuffers) {
        SGEGSpectralGroup sg;

        sg.root_name.resize(fb->root_name.size() + 1);
        std::memcpy(sg.root_name.data(), fb->root_name.c_str(), sg.root_name.size() * sizeof(char));
        sg.wavelengths = fb->wavelengths_nm;

        std::pair<int, int> n_bits_dc;

        quantization_from_exr(fb->pixel_type, n_bits_dc);

        std::vector<std::vector<float>> compressed_moments;
        std::vector<uint8_t> relative_scales;
        std::vector<std::pair<int, int>> quantization_curve;
        std::vector<uint32_t> downsampling_factor_curve;
        std::vector<float> compression_curve;

        stats_data quantization_error, compression_error;

        compress_spectral_framebuffer(
            method,
            fb,
            exr_in.width(), exr_in.height(),
            n_bits_dc,      n_bits_ac1,      use_flat_quantization, quantization_curve,
            normalize_moments,
            compression_dc, compression_ac1, use_flat_compression,  compression_curve, compression_effort,
            downsampling_factor_ac1, downsampling_factor_curve,
            compressed_moments,
            sg,
            log_is_active,
            log_content,
            quantization_error,
            compression_error
        );

        quantization_errors.push_back(quantization_error);
        compression_errors.push_back(compression_error);

        assert(quantization_curve.size() == compressed_moments.size());
        assert(compression_curve.size() == compressed_moments.size());

        // Now we can save to JPEG XL
        for (size_t m = 0; m < compressed_moments.size(); m++) {
            const size_t idx = jxl_out.appendFramebuffer(
                compressed_moments[m],
                1,
                quantization_curve[m],
                downsampling_factor_curve[m],
                compression_curve[m],
                fb->root_name.c_str());

            sg.layer_indices.push_back(idx);
        }

        box.spectral_groups.push_back(sg);
    }

    // Append the extra framebuffers as it
    for (const GreyFramebuffer* fb: extra_framebuffers) {
        SGEGGrayGroup gg;

        gg.layer_name.resize(fb->layer_name.size() + 1);
        std::memcpy(gg.layer_name.data(), fb->layer_name.c_str(), gg.layer_name.size() * sizeof(char));

        std::pair<int, int> n_bits;

        quantization_from_exr(fb->pixel_type, n_bits);

        gg.layer_index = jxl_out.appendFramebuffer(
            fb->image_data,
            1,
            n_bits,
            1,
            compression_dc,
            fb->layer_name.c_str()
        );

        box.gray_groups.push_back(gg);
    }

    jxl_out.setBox(box);
    jxl_out.write(filename_out, compression_effort);

    if (has_dump_file) {
        jxl_out.dump(filename_dump);
    }

    clock_end = std::chrono::steady_clock::now();

    // Logging stuff

    if (log_is_active) {
        auto diff = clock_end - clock_start;
        log_content << "Total duration: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

        std::ofstream log_file(log_filepath);
        log_file << log_content.str();
    }

    // Saves log in binary form
    if (has_bin_log) {
        FILE *f = fopen(binarylog_filepath.c_str(), "wb");

        uint32_t n = quantization_errors.size();

        fwrite(&n, sizeof(uint32_t), 1, f);

        for (size_t i = 0; i < quantization_errors.size(); i++) {
            fwrite(&quantization_errors[i].rmse_error, sizeof(double), 1, f);
            fwrite(&quantization_errors[i].rrmse_error, sizeof(double), 1, f);
            fwrite(&quantization_errors[i].max_error, sizeof(double), 1, f);

            fwrite(&compression_errors[i].rmse_error, sizeof(double), 1, f);
            fwrite(&compression_errors[i].rrmse_error, sizeof(double), 1, f);
            fwrite(&compression_errors[i].max_error, sizeof(double), 1, f);
        }

        fclose(f);
    }

    return 0;
}

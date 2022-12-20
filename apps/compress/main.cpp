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
#include <sstream>
#include <limits>
#include <regex>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <fstream>
#include <cassert>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include <tclap/CmdLine.h>

#include <JXLImage.h>
#include <EXRSpectralImage.h>

#include <moments.h>
#include <moments_image.h>
#include <curve_quantization.h>
#include <curve_compression.h>

/**
 * TODO:
 * - Add an argument
 *    - to control the number of bits for starting the quantization curves
 *    - to control the JXL quality factor -> TODO, the option is there but
 *      seems not affecting file size while decreasing quality. Also 1 is
 *      supposed to look the same but does drastically decrease the quality
 *    - to control spatial downsampling (needs as well serious work in the JXL class)
 * - Fix runtime error when running with CAVE database
 * - RGB layers are only computed on the root (maybe a desirable behaviour though...)
 */


void compress_spectral_framebuffer(
    SpectralStorageMethod method,
    const SpectralFramebuffer* framebuffer,
    uint32_t width, uint32_t height,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve,
    float start_compression_curve,
    std::vector<float>& compression_curve,
    std::vector<std::vector<float>>& compressed_moments,
    SGEGSpectralGroup& sg,
    bool log,
    std::stringstream& log_stream)
{
    sg.method = method;
    double quantization_error;
    double compression_error;

    std::chrono::time_point<std::chrono::steady_clock>  clock_start, clock_end;

    std::vector<double> compressed_moments_d;
    std::vector<uint8_t> relative_scales;
    std::vector<double> mins_d, maxs_d;
    double global_min_d, global_max_d;

    const uint32_t n_moments = framebuffer->wavelengths_nm.size();
    const uint32_t n_pixels = framebuffer->image_data.size() / framebuffer->wavelengths_nm.size();

    // TODO: this shall be revisited
    // TODO: this is invalid for polarization layers that can be negative
    // Condition buffers to ensure everything is in order
    std::vector<double> spectral_wavelengths(framebuffer->wavelengths_nm.size());
    std::vector<double> spectral_framebuffer(framebuffer->image_data.size());

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_wavelengths.size(); i++) {
        spectral_wavelengths[i] = framebuffer->wavelengths_nm[i];
    }

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_framebuffer.size(); i++) {
        double v = framebuffer->image_data[i];

        if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
            v = 1e-8;
        }

        spectral_framebuffer[i] = v;
    }

    // TODO: this has to be removed and calculated on a per compression scheme basis
    compression_curve.resize(n_moments);

    for (size_t m = 0; m < n_moments; m++) {
        compression_curve[m] = start_compression_curve;
    }

    std::cout << "[WARNING] The compression curve may be not computed!" << std::endl;

    switch (sg.method) {
        case LINEAR:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = linear_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc, n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            compression_error = linear_compute_compression_curve(
                spectral_wavelengths, spectral_framebuffer,
                width, height, n_moments,
                quantization_curve,
                start_compression_curve,
                compression_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            linear_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            if (log) {
                log_stream << "Compression curve:" << std::endl;

                for (const float& c : compression_curve) {
                    log_stream << c << " ";
                }

                log_stream << std::endl << compression_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }
            break;

        case BOUNDED:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = bounded_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc, n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            bounded_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }
            break;

        case UNBOUNDED:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = unbounded_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc, n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            unbounded_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }
            break;

        case UNBOUNDED_TO_BOUNDED:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = unbounded_to_bounded_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc, n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            unbounded_to_bounded_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }
            break;

        case UPPERBOUND:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = upperbound_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc,
                n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            upperbound_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d,
                relative_scales,
                global_max_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments + 1);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            // Copy the local scaling
            quantization_curve.push_back(8);
            compressed_moments[n_moments].resize(n_pixels);

            for (size_t px = 0; px < n_pixels; px++) {
                compressed_moments[n_moments][px] = (double)relative_scales[px] / (double)std::numeric_limits<uint8_t>::max();
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }

            sg.global_max = global_max_d;
            break;

        case TWOBOUNDS:
            if (log) {
                clock_start = std::chrono::steady_clock::now();
            }

            // Create a quantization profil
            quantization_error = twobounds_compute_quantization_curve(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                n_bits_dc,
                n_bits_ac1,
                quantization_curve
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Quantization: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;

                clock_start = std::chrono::steady_clock::now();
            }

            twobounds_compress_spectral_image(
                spectral_wavelengths, spectral_framebuffer,
                n_pixels, n_moments,
                compressed_moments_d,
                mins_d, maxs_d,
                relative_scales,
                global_min_d,
                global_max_d
            );

            if (log) {
                clock_end = std::chrono::steady_clock::now();
                auto diff = clock_end - clock_start;
                log_stream << "Compression: ";
                log_stream << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;;
            }

            if (log) {
                log_stream << "Quantization curve:" << std::endl;

                for (const int& q : quantization_curve) {
                    log_stream << q << " ";
                }

                log_stream << std::endl << quantization_error << std::endl;
            }

            // Copy back and implicit conversion to float
            compressed_moments.resize(n_moments + 1);

            for (size_t m = 0; m < n_moments; m++) {
                compressed_moments[m].resize(n_pixels);

                for (size_t px = 0; px < n_pixels; px++) {
                    compressed_moments[m][px] = compressed_moments_d[n_moments * px + m];
                }
            }

            // Copy the local scaling
            quantization_curve.push_back(8);
            compressed_moments[n_moments].resize(n_pixels);

            for (size_t px = 0; px < n_pixels; px++) {
                compressed_moments[n_moments][px] = (double)relative_scales[px] / (double)std::numeric_limits<uint8_t>::max();
            }

            sg.mins.resize(n_moments - 1);
            sg.maxs.resize(n_moments - 1);

            for (size_t m = 0; m < n_moments - 1; m++) {
                sg.mins[m] = mins_d[m];
                sg.maxs[m] = maxs_d[m];
            }

            sg.global_min = global_min_d;
            sg.global_max = global_max_d;
            break;
    }
}


void quantization_from_exr(PixelType type, size_t& n_bits, size_t& n_exponent_bits) {
    switch (type) {
        case PixelType::UINT:
            n_bits = 32;
            n_exponent_bits = 0;
            break;

        case PixelType::HALF:
            n_bits = 16;
            n_exponent_bits = 5;
            break;

        case PixelType::FLOAT:
            n_bits = 32;
            n_exponent_bits = 8;
            break;

        default:
            throw std::runtime_error("Unknown pixel type");
            break;
    }
}


class FrameDistanceConstraint: public TCLAP::Constraint<float>
{
public:
    virtual std::string description() const
    {
        return "Sets the distance level for lossy compression";
    }

    virtual std::string shortID() const
    {
        return "FrameDistanceConstraint";
    }

    virtual bool check(const float &value) const
    {
        return (value >= 0.f) && (value <= 15.f);
    }
};


int main(int argc, char *argv[])
{
    std::string filename_in, filename_out;
    float frame_distance = .1f;
    int n_bits_start_quantization = 10;
    SpectralStorageMethod method = TWOBOUNDS;
    std::chrono::time_point<std::chrono::steady_clock>  clock_start, clock_end;
    bool log_is_active = false;
    std::string log_filepath;
    std::stringstream log_content;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compress a spectral image");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg("Input", "Spectral EXR input file", true, "input.exr", "path");
        cmd.add(inputFileArg);

        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "JPEG XL output file", true, "output.jxl", "path");
        cmd.add(outputFileArg);

        TCLAP::ValueArg<std::string> logFileArg("l", "log", "File to log into", false, "log.txt", "path");
        cmd.add(logFileArg);

        FrameDistanceConstraint frameDistanceConstraint;
        TCLAP::ValueArg<float> frameDistanceArg("d", "frame_distance", "Distance level for lossy compression (compression rate). Range: 0 .. 15", false, .1f, &frameDistanceConstraint);
        cmd.add(frameDistanceArg);

        TCLAP::ValueArg<int> quantizationStartArg("q", "quantization", "Starting number of bits for quantizing the first AC component.", false, 10, "integer");
        cmd.add(quantizationStartArg);

        std::vector<std::string> allowedCompressionMethods;
        allowedCompressionMethods.push_back("linear");
        allowedCompressionMethods.push_back("bounded");
        allowedCompressionMethods.push_back("unbounded");
        allowedCompressionMethods.push_back("unbounded_to_bounded");
        allowedCompressionMethods.push_back("upperbound");
        allowedCompressionMethods.push_back("twobounds");
        TCLAP::ValuesConstraint<std::string> allowedCompressionMethodsVals(allowedCompressionMethods);
        TCLAP::ValueArg<std::string> momentCompressionMethodArg("m", "method", "Moment compression method to use", false, "twobounds", &allowedCompressionMethodsVals);
        cmd.add(momentCompressionMethodArg);

        cmd.parse(argc, argv);

        filename_in    = inputFileArg.getValue();
        filename_out   = outputFileArg.getValue();
        frame_distance = frameDistanceArg.getValue();
        n_bits_start_quantization = quantizationStartArg.getValue();

        log_is_active = logFileArg.isSet();

        if (log_is_active) {
            log_filepath = logFileArg.getValue();
        }

        if (momentCompressionMethodArg.getValue() == "linear") {
            method = LINEAR;
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
    }

    // Timing setup
    clock_start = std::chrono::steady_clock::now();

    // Load input file
    EXRSpectralImage exr_in(filename_in);

    const std::vector<SpectralFramebuffer*>& spectral_framebuffers = exr_in.getSpectralFramebuffers();
    const std::vector<GreyFramebuffer*>& extra_framebuffers = exr_in.getExtraFramebuffers();

    // Create output file
    JXLImage jxl_out(exr_in.width(), exr_in.height());
    SGEGBox box;

    box.exr_attributes = exr_in.getAttributesData();

    // Run compression for each spectral group
    for (const SpectralFramebuffer* fb: spectral_framebuffers) {
        SGEGSpectralGroup sg;

        sg.root_name.resize(fb->root_name.size() + 1);
        std::memcpy(sg.root_name.data(), fb->root_name.c_str(), sg.root_name.size() * sizeof(char));
        sg.wavelengths = fb->wavelengths_nm;

        size_t main_n_bits;
        size_t main_n_exponent_bits;

        quantization_from_exr(fb->pixel_type, main_n_bits, main_n_exponent_bits);

        std::vector<std::vector<float>> compressed_moments;
        std::vector<uint8_t> relative_scales;
        std::vector<int> quantization_curve;
        std::vector<float> compression_curve;

        compress_spectral_framebuffer(
            method,
            fb,
            exr_in.width(), exr_in.height(),
            main_n_bits,
            n_bits_start_quantization,
            quantization_curve,
            frame_distance,
            compression_curve,
            compressed_moments,
            sg,
            log_is_active,
            log_content
        );

        assert(quantization_curve.size() == compressed_moments.size());
        assert(compression_curve.size() == compressed_moments.size());

        // Now we can save to JPEG XL
        for (size_t m = 0; m < compressed_moments.size(); m++) {
            float n_bits;
            float n_exponent_bits;

            if (m == 0) {
                n_bits          = main_n_bits;
                n_exponent_bits = main_n_exponent_bits;
            } else {
                n_bits          = quantization_curve[m];
                n_exponent_bits = 0;
            }

            const size_t idx = jxl_out.appendFramebuffer(
                compressed_moments[m],
                1,
                n_bits,
                n_exponent_bits,
                1,
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

        size_t n_bits;
        size_t n_exponent_bits;

        quantization_from_exr(fb->pixel_type, n_bits, n_exponent_bits);

        gg.layer_index = jxl_out.appendFramebuffer(
            fb->image_data,
            1,
            n_bits,
            n_exponent_bits,
            1,
            frame_distance,
            fb->layer_name.c_str()
        );

        box.gray_groups.push_back(gg);
    }

    jxl_out.setBox(box);
    jxl_out.write(filename_out);

    clock_end = std::chrono::steady_clock::now();

    if (log_is_active) {
        auto diff = clock_end - clock_start;
        log_content << "Total duration: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

        std::ofstream log_file(log_filepath);
        log_file << log_content.str();
    }

#ifndef NDEBUG
    // Test dump
    jxl_out.dump("jxl_dump");
    exr_in.dump("exr_dump");
#endif // NDEBUG

    return 0;
}

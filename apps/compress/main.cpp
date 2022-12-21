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
 *    - to control spatial downsampling (needs as well serious work in the JXL class)
 * - RGB layers are only computed on the root (maybe a desirable behaviour though...)
 */


void generate_quantization_curve(
    SpectralStorageMethod method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    bool uses_constant_quantization,
    std::vector<int>& quantization_curve,
    double& timing)
{
    if (uses_constant_quantization) {
        quantization_curve.resize(n_moments);
        quantization_curve[0] = n_bits_dc;

        for (size_t i = 1; i < n_moments; i++) {
            quantization_curve[i] = n_bits_ac1;
        }
        timing = 0;
    } else {
        auto clock_start = std::chrono::steady_clock::now();

        switch (method) {
            case LINEAR:
                linear_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case BOUNDED:
                bounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UNBOUNDED:
                unbounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UNBOUNDED_TO_BOUNDED:
                unbounded_to_bounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UPPERBOUND:
                upperbound_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case TWOBOUNDS:
                twobounds_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
        }

        auto clock_end = std::chrono::steady_clock::now();
        timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();
    }

    assert(quantization_curve.size() == n_moments);
}


void generate_compression_curve(
    SpectralStorageMethod method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    bool uses_constant_compression,
    std::vector<float>& compression_curve,
    double& timing)
{
    if (uses_constant_compression) {
        compression_curve.resize(n_moments);
        compression_curve[0] = compression_dc;

        for (size_t i = 1; i < n_moments; i++) {
            compression_curve[i] = compression_ac1;
        }

        timing = 0;
        // TODO error
    } else {
        auto clock_start = std::chrono::steady_clock::now();

        switch (method) {
            case LINEAR:
                linear_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
            case BOUNDED:
                bounded_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
            case UNBOUNDED:
                unbounded_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
            case UNBOUNDED_TO_BOUNDED:
                unbounded_to_bounded_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
            case UPPERBOUND:
                upperbound_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
            case TWOBOUNDS:
                twobounds_compute_compression_curve(
                    wavelengths, spectral_image,
                    width, height, n_moments,
                    quantization_curve,
                    compression_dc, compression_ac1,
                    compression_curve
                );
                break;
        }

        auto clock_end = std::chrono::steady_clock::now();
        timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();
    }

    assert(compression_curve.size() == n_moments);
}


template<typename T>
void compress_image(
    SpectralStorageMethod method,
    const std::vector<double>& wavelengths,
    const std::vector<double> spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    std::vector<double>& compressed_moments,
    std::vector<double>& mins, std::vector<double>& maxs,
    std::vector<T>& relative_scales,
    double& global_min, double& global_max,
    double& timing)
{
    const size_t n_pixels = width * height;

    auto clock_start = std::chrono::steady_clock::now();

    switch (method) {
        case LINEAR:
            linear_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs
            );
            break;

        case BOUNDED:
            bounded_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs
            );
            break;

        case UNBOUNDED:
            unbounded_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs
            );
            break;

        case UNBOUNDED_TO_BOUNDED:
            unbounded_to_bounded_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs
            );
            break;

        case UPPERBOUND:
            upperbound_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs,
                relative_scales,
                global_max
            );
            break;

        case TWOBOUNDS:
            twobounds_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs,
                relative_scales,
                global_min,
                global_max
            );
            break;
    }

    auto clock_end = std::chrono::steady_clock::now();
    timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();

    assert(compressed_moments_d.size() == n_pixels * n_moments);
    assert(mins_d.size() == n_moments - 1);
    assert(maxs_d.size() == n_moments - 1);

#ifndef NDEBUG
    if (method == UPPERBOUND || method == TWOBOUNDS) {
        assert(relative_scales.size() == n_pixels);
    }
#endif // NDEBUG
}


double error_for_quantization_curve(
    SpectralStorageMethod method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max)
{
    switch(method) {
        case LINEAR:
            return linear_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case BOUNDED:
            return bounded_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UNBOUNDED:
            return unbounded_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UNBOUNDED_TO_BOUNDED:
            return unbounded_to_bounded_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UPPERBOUND:
            return upperbound_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_max,
                quantization_curve
            );
        case TWOBOUNDS:
            return twobounds_error_for_quantization_curve(
                wavelengths, spectral_image,
                width * height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_min, global_max,
                quantization_curve
            );
    }

    assert(0);
    return 0;
}


double error_for_quantization_and_compression_curves(
    SpectralStorageMethod method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max)
{
    switch(method) {
        case LINEAR:
            return linear_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve
            );
        case BOUNDED:
            return bounded_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve
            );
        case UNBOUNDED:
            return unbounded_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve
            );
        case UNBOUNDED_TO_BOUNDED:
            return unbounded_to_bounded_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve
            );
        case UPPERBOUND:
            return upperbound_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_max,
                quantization_curve,
                compression_curve
            );
        case TWOBOUNDS:
            return twobounds_error_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_min, global_max,
                quantization_curve,
                compression_curve
            );
    }

    assert(0);
    return 0;
}


void compress_spectral_framebuffer(
    SpectralStorageMethod method,
    const SpectralFramebuffer* framebuffer,
    uint32_t width, uint32_t height,
    int n_bits_dc,
    int n_bits_ac1,
    bool uses_constant_quantization,
    std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    bool uses_constant_compression,
    std::vector<float>& compression_curve,
    std::vector<std::vector<float>>& compressed_moments,
    SGEGSpectralGroup& sg,
    bool log,
    std::stringstream& log_stream)
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
    std::vector<double> spectral_image(framebuffer->image_data.size());

    #pragma omp parallel for
    for (size_t i = 0; i < wavelengths.size(); i++) {
        wavelengths[i] = framebuffer->wavelengths_nm[i];
    }

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_image.size(); i++) {
        double v = framebuffer->image_data[i];

        if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
            v = 1e-8;
        }

        spectral_image[i] = v;
    }

    double quantization_curve_timing(0),
           compression_curve_timing(0),
           compression_timing(0);

    double quantization_error(0),
           compression_error(0);

    generate_quantization_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        n_bits_dc, n_bits_ac1,
        uses_constant_quantization,
        quantization_curve,
        quantization_curve_timing
    );

    generate_compression_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        compression_dc, compression_ac1,
        uses_constant_compression,
        compression_curve,
        compression_curve_timing
    );

    compress_image(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        compressed_moments_d,
        mins_d, maxs_d,
        relative_scales,
        global_min_d, global_max_d,
        compression_timing
    );

    // Compute errors (may be done twice if the quantization / compression
    // profile is built at runtime)

    quantization_error = error_for_quantization_curve(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        compressed_moments_d, mins_d, maxs_d,
        relative_scales, global_min_d, global_max_d
    );

    compression_error = error_for_quantization_and_compression_curves(
        method,
        wavelengths, spectral_image,
        width, height, n_moments,
        quantization_curve,
        compression_curve,
        compressed_moments_d, mins_d, maxs_d,
        relative_scales, global_min_d, global_max_d
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

        quantization_curve.push_back(8);
        quantization_curve.push_back(compression_dc);

        compressed_moments.resize(n_moments + 1);
        compressed_moments[n_moments].resize(n_pixels);

        for (size_t px = 0; px < n_pixels; px++) {
            compressed_moments[n_moments][px] = (float)((double)relative_scales[px] / (double)std::numeric_limits<uint8_t>::max());
        }

        sg.global_min = global_min_d;
        sg.global_max = global_max_d;
    }

    if (log) {
        log_stream << "Quantization curve:" << std::endl;
        for (const int& q : quantization_curve) {
            log_stream << q << " ";
        }
        log_stream << std::endl << quantization_error << std::endl;

        log_stream << "Compression curve:" << std::endl;
        for (const float& c : compression_curve) {
            log_stream << c << " ";
        }
        log_stream << std::endl << compression_error << std::endl;

        log_stream << "Timings:" << std::endl;
        log_stream << "Quantization curve: " << quantization_curve_timing << " ms" << std::endl;
        log_stream << "Compression curve:  " << compression_curve_timing << " ms" << std::endl;
        log_stream << "Compression:        " << compression_timing << " ms" << std::endl;
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
        return "0..15";
    }

    virtual bool check(const float &value) const
    {
        return (value >= 0.f) && (value <= 15.f);
    }
};


int main(int argc, char *argv[])
{
    std::string filename_in, filename_out;
    bool log_is_active = false;
    std::string log_filepath;
    std::stringstream log_content;

    float compression_dc = .1f;
    float compression_ac1 = .1f;
    bool use_flat_compression = false;

    int n_bits_ac1 = 10;
    bool use_flat_quantization = false;

    SpectralStorageMethod method = TWOBOUNDS;

    std::chrono::time_point<std::chrono::steady_clock>  clock_start, clock_end;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Compress a spectral image");

        // Input / output
        TCLAP::UnlabeledValueArg<std::string> inputFileArg("Input", "Spectral EXR input file", true, "input.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "JPEG XL output file", true, "output.jxl", "path");
        TCLAP::ValueArg<std::string> logFileArg("l", "log", "File to log into", false, "log.txt", "path");
        cmd.add(inputFileArg);
        cmd.add(outputFileArg);
        cmd.add(logFileArg);

        // Compresion tweaking
        FrameDistanceConstraint frameDistanceConstraint;
        TCLAP::ValueArg<float> frameDistanceDCArg("a", "frame_distance_dc", "Distance level for lossy compression (compression rate) on the DC component.", false, .1f, &frameDistanceConstraint);
        TCLAP::ValueArg<float> frameDistanceACArg("b", "frame_distance_ac", "Distance level for lossy compression (compression rate) on the first AC component.", false, .1f, &frameDistanceConstraint);
        TCLAP::SwitchArg useFlatCompressionArg("c", "c_flat", "Use a flat compression curve");
        cmd.add(frameDistanceDCArg);
        cmd.add(frameDistanceACArg);
        cmd.add(useFlatCompressionArg);

        // Quantization tweaking
        TCLAP::ValueArg<int> quantizationStartArg("q", "quantization", "Starting number of bits for quantizing the first AC component.", false, 10, "integer");
        TCLAP::SwitchArg useFlatQuantizationArg("u", "q_flat", "Use a flat quantization curve");
        cmd.add(quantizationStartArg);
        cmd.add(useFlatQuantizationArg);

        // Moment storage method tweaking
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

        filename_in     = inputFileArg.getValue();
        filename_out    = outputFileArg.getValue();
        log_is_active   = logFileArg.isSet();

        if (log_is_active) {
            log_filepath = logFileArg.getValue();
        }

        compression_dc       = frameDistanceDCArg.getValue();
        compression_ac1      = frameDistanceACArg.getValue();
        use_flat_compression = useFlatCompressionArg.getValue();

        n_bits_ac1            = quantizationStartArg.getValue();
        use_flat_quantization = useFlatQuantizationArg.getValue();

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

        size_t n_bits_dc;
        size_t main_n_exponent_bits;

        quantization_from_exr(fb->pixel_type, n_bits_dc, main_n_exponent_bits);

        std::vector<std::vector<float>> compressed_moments;
        std::vector<uint8_t> relative_scales;
        std::vector<int> quantization_curve;
        std::vector<float> compression_curve;

        compress_spectral_framebuffer(
            method,
            fb,
            exr_in.width(), exr_in.height(),
            n_bits_dc,      n_bits_ac1,      use_flat_quantization, quantization_curve,
            compression_dc, compression_ac1, use_flat_compression,  compression_curve,
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
                n_bits          = n_bits_dc;
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
            compression_dc,
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

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

#include "curve_compression.h"
#include "curve_quantization.h"

#include "moments_image.h"
#include "moments_error.h"
#include "Util.h"

#include <jxl/encode.h>
#include <jxl/encode_cxx.h>
#include <jxl/decode.h>
#include <jxl/decode_cxx.h>
#include <jxl/thread_parallel_runner.h>
#include <jxl/thread_parallel_runner_cxx.h>

#include <half.h>

#include <vector>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <chrono>

// #define DEBLOG

// ----------------------------------------------------------------------------

#define CHECK_JXL_ENC_STATUS(status)                                          \
    if (JXL_ENC_SUCCESS != (status)) {                                        \
        std::stringstream err_msg;                                            \
        err_msg << "[ERR] JXL_ENCODER: "                                      \
                << __FILE__ << ":" << __LINE__                                \
                << "> " << status;                                            \
        throw std::runtime_error(err_msg.str());                              \
    }                                                                         \


#define CHECK_JXL_DEC_STATUS(status)                                          \
    if (JXL_DEC_SUCCESS != (status)) {                                        \
        std::stringstream err_msg;                                            \
        err_msg << "[ERR] JXL_DECODER: "                                      \
                << __FILE__ << ":" << __LINE__                                \
                << "> " << status;                                            \
        throw std::runtime_error(err_msg.str());                              \
    }                                                                         \

// ----------------------------------------------------------------------------

void compress_decompress_framebuffer(
    const std::vector<float>& framebuffer_in,
    std::vector<float>& framebuffer_out,
    uint32_t width, uint32_t height,
    std::pair<int, int> bits_per_sample,
    uint32_t subsampling_factor,
    float frame_distance,
    int effort)
{
    assert(framebuffer_in.size() == width * height);
    framebuffer_out.resize(framebuffer_in.size());

    // 1. Compress the framebuffer using JXL
    std::vector<float> quantized_framebuffer;

    JxlEncoderStatus           enc_status;
    JxlEncoderPtr              enc;
    JxlThreadParallelRunnerPtr runner;
    JxlBasicInfo               basic_info;
    JxlColorEncoding           color_encoding;

    runner = JxlThreadParallelRunnerMake(
            nullptr,
            JxlThreadParallelRunnerDefaultNumWorkerThreads()
        );

    enc = JxlEncoderMake(nullptr);

    enc_status = JxlEncoderSetParallelRunner(
        enc.get(),
        JxlThreadParallelRunner,
        runner.get()
    );
    CHECK_JXL_ENC_STATUS(enc_status);

    JxlEncoderInitBasicInfo(&basic_info);

    basic_info.xsize                    = width;
    basic_info.ysize                    = height;
    basic_info.num_extra_channels       = 0;
    basic_info.num_color_channels       = 1;
    basic_info.bits_per_sample          = bits_per_sample.first;
    basic_info.exponent_bits_per_sample = bits_per_sample.second;
    basic_info.uses_original_profile    = JXL_TRUE;

    enc_status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);
    CHECK_JXL_ENC_STATUS(enc_status);

    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, effort);
    CHECK_JXL_ENC_STATUS(enc_status);

    const bool encodes_lossless = frame_distance == 0;

    if (encodes_lossless) {
        enc_status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    } else {
        enc_status = JxlEncoderSetFrameLossless(frame_settings, JXL_FALSE);
        enc_status = JxlEncoderSetFrameDistance(frame_settings, frame_distance);
    }
    CHECK_JXL_ENC_STATUS(enc_status);

    enc_status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, subsampling_factor);
    CHECK_JXL_ENC_STATUS(enc_status);

    color_encoding.color_space       = JXL_COLOR_SPACE_GRAY;
    color_encoding.white_point       = JXL_WHITE_POINT_D65;
    color_encoding.primaries         = JXL_PRIMARIES_SRGB;
    color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
    color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

    enc_status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);
    CHECK_JXL_ENC_STATUS(enc_status);

    JxlPixelFormat format;
    format.num_channels = 1;
    format.data_type    = JXL_TYPE_FLOAT,
    format.endianness   = JXL_NATIVE_ENDIAN;
    format.align        = 0;

    // FIXME: There is a bug in the current JXL implementation: it ignores
    //        the quantization when the compression is not lossless.
    //        This issue is addressed with a hack.
    if (encodes_lossless || basic_info.exponent_bits_per_sample != 0) {
        // When we are dealing with lossless file and halfs, libjxl
        // triggers an error if the half rounding cause imprecision.
        // So we explicitely cast to half and back to float before
        // providing the buffer.
        if (basic_info.exponent_bits_per_sample == 5 && basic_info.bits_per_sample == 16) {
            quantized_framebuffer.resize(width * height);

            #pragma omp parallel for
            for (size_t i = 0; i < width * height; i++) {
                quantized_framebuffer[i] =
                    imath_half_to_float(
                        imath_float_to_half(framebuffer_in[i]
                    )
                );
            }

            enc_status = JxlEncoderAddImageFrame(
                frame_settings,
                &format,
                quantized_framebuffer.data(),
                width * height * sizeof(float)
            );
        } else {
            enc_status = JxlEncoderAddImageFrame(
                frame_settings,
                &format,
                framebuffer_in.data(),
                width * height * sizeof(float)
            );
        }
    } else {
        quantized_framebuffer.resize(width * height);

        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            quantized_framebuffer[i] = Util::quantize_dequantize(
                framebuffer_in[i],
                basic_info.bits_per_sample
            );
        }

        enc_status = JxlEncoderAddImageFrame(
            frame_settings,
            &format,
            quantized_framebuffer.data(),
            width * height * sizeof(float)
        );
    }

    CHECK_JXL_ENC_STATUS(enc_status);

    JxlEncoderCloseInput(enc.get());

    std::vector<uint8_t> compressed(64);

    uint8_t* next_out = compressed.data();
    size_t avail_out  = compressed.size();

    enc_status = JXL_ENC_NEED_MORE_OUTPUT;

    while (enc_status == JXL_ENC_NEED_MORE_OUTPUT) {
        enc_status = JxlEncoderProcessOutput(enc.get(), &next_out, &avail_out);

        if (enc_status == JXL_ENC_NEED_MORE_OUTPUT) {
            size_t offset = next_out - compressed.data();
            compressed.resize(compressed.size() * 2);

            next_out  = compressed.data() + offset;
            avail_out = compressed.size() - offset;
        }
    }

    CHECK_JXL_ENC_STATUS(enc_status);

    compressed.resize(next_out - compressed.data());

    // Now we have the compressed data in `compresed`

    // 2. Decompress the framebuffer using JXL
    JxlDecoderStatus dec_status;
    JxlDecoderPtr    dec;

    dec = JxlDecoderMake(nullptr);
    dec_status = JxlDecoderSetParallelRunner(
        dec.get(),
        JxlThreadParallelRunner,
        runner.get()
    );
    CHECK_JXL_DEC_STATUS(dec_status);

    dec_status = JxlDecoderSubscribeEvents(
        dec.get(),
        JXL_DEC_BASIC_INFO |
        JXL_DEC_FULL_IMAGE);
    CHECK_JXL_DEC_STATUS(dec_status);

    JxlDecoderSetInput(dec.get(), compressed.data(), compressed.size());

    for (JxlDecoderStatus status_process = JXL_DEC_NEED_MORE_INPUT;
        status_process != JXL_DEC_FULL_IMAGE;
        status_process = JxlDecoderProcessInput(dec.get()))
    {
        switch (status_process) {
            case JXL_DEC_BASIC_INFO:
                {
                    // Read metadata and allocate memory
                    JxlBasicInfo basic_info;
                    // JxlExtraChannelInfo extra_info;
                    std::vector<char> layer_name;

                    // Main layer metadata
                    dec_status = JxlDecoderGetBasicInfo(dec.get(), &basic_info);

                    // We do not expect any extra framebuffer
                    assert(basic_info.num_color_channels == 1);
                    assert(basic_info.xsize == width);
                    assert(basic_info.ysize == height);
                }
                break;

            case JXL_DEC_NEED_IMAGE_OUT_BUFFER:
                {
                    JxlPixelFormat format;
                    format.num_channels = 1;
                    format.data_type    = JXL_TYPE_FLOAT,
                    format.endianness   = JXL_NATIVE_ENDIAN;
                    format.align        = 0;

                    dec_status = JxlDecoderSetImageOutBuffer(
                        dec.get(),
                        &format,
                        framebuffer_out.data(),
                        width * height * sizeof(float)
                    );
                    CHECK_JXL_DEC_STATUS(dec_status);
                }
                break;
            case JXL_DEC_NEED_MORE_INPUT:
                break;
            case JXL_DEC_ERROR:
                // CHECK_JXL_DEC_STATUS(status_process);
                // break;
            default:
                std::stringstream err_msg;
                err_msg << "Unknown decoder status: " << status_process;
                throw std::runtime_error(err_msg.str());
                break;
        }
    }

    JxlDecoderReleaseInput(dec.get());
    CHECK_JXL_DEC_STATUS(dec_status);
}


void compress_decompress_single_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    std::pair<int, int> bits_per_sample,
    uint32_t subsampling_factor,
    size_t i, float frame_distance,
    int effort)
{
    assert(input_image.size() == width * height * n_moments);

    // Extract a single moment from the input image
    std::vector<float> input_framebuffer(width * height);
    std::vector<float> compressed_framebuffer;

    #pragma omp parallel for
    for (size_t px = 0; px < width * height; px++) {
        input_framebuffer[px] = (float)input_image[px * n_moments + i];
    }

    compress_decompress_framebuffer(
        input_framebuffer,
        compressed_framebuffer,
        width, height,
        bits_per_sample,
        subsampling_factor,
        frame_distance,
        effort
    );

    // Copy back the data to the output image
    output_image.resize(width * height * n_moments);

    std::memcpy(
        output_image.data(),
        input_image.data(),
        width * height * n_moments * sizeof(double)
    );

    // Change what has changed
    #pragma omp parallel for
    for (size_t px = 0; px < width * height; px++) {
        output_image[px * n_moments + i] = (double)compressed_framebuffer[px];
    }
}


void compress_decompress_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(input_image.size() == width * height * n_moments);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    output_image.resize(width * height * n_moments);

    // JXL compression for every layers
    #pragma omp parallel for
    for (size_t m = 0; m < n_moments; m++) {
        std::vector<float> framebuffer_in(width * height);
        std::vector<float> framebuffer_out;

        // Copy cast to float
        #pragma omp parallel for
        for (size_t px = 0; px < width * height; px++) {
            framebuffer_in[px] = input_image[px * n_moments + m];
        }

        // Compress / decompress
        compress_decompress_framebuffer(
            framebuffer_in,
            framebuffer_out,
            width, height,
            quantization_curve[m],
            subsampling_factor_curve[m],
            compression_curve[m],
            effort);

        // Copy cast to float
        #pragma omp parallel for
        for (size_t px = 0; px < width * height; px++) {
            output_image[px * n_moments + m] = framebuffer_out[px];
        }
    }
}


/*****************************************************************************/
/* Create compression curves                                                 */
/*****************************************************************************/

float sigmoid(float x)
{
    return 1.f / (1.f + std::exp(-x));
}


float deterministic_compression(float fd_start, float fd_end, size_t n_moments, float x)
{
    const float n = (float) n_moments;
    const float s_0 = sigmoid(10.f * (2.f - n)/n - 1.f);

    return fd_start
        + (fd_end - fd_start) *
             (sigmoid(10.f * (2.f *     x     - n) / n - 1.f) - s_0)
           / (sigmoid(10.f * (2.f * (n - 1.f) - n) / n - 1.f) - s_0);
}


void compute_compression_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    CompressionCurveType compression_curve_type,
    std::vector<float>& compression_curve,
    float effort,
    double& timing)
{
    switch (compression_curve_type) {
        case COMPRESSION_FLAT:
            compression_curve.resize(n_moments);
            compression_curve[0] = compression_dc;

            for (size_t i = 1; i < n_moments; i++) {
                compression_curve[i] = compression_ac1;
            }

            timing = 0;
            break;
        case COMPRESSION_DYNAMIC:
            {
                auto clock_start = std::chrono::steady_clock::now();

                switch (method) {
                    case LINEAR:
                        linear_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case LINAVG:
                        linavg_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case BOUNDED:
                        bounded_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case UNBOUNDED:
                        unbounded_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case UNBOUNDED_TO_BOUNDED:
                        unbounded_to_bounded_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case UPPERBOUND:
                        upperbound_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                    case TWOBOUNDS:
                        twobounds_compute_compression_curve(
                            wavelengths, spectral_image,
                            width, height, n_moments,
                            quantization_curve,
                            normalize_moments,
                            subsampling_factor_curve,
                            compression_dc, compression_ac1,
                            compression_curve,
                            effort
                        );
                        break;
                }

                auto clock_end = std::chrono::steady_clock::now();
                timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();
            }
            break;

        case COMPRESSION_DETERMINISTIC:
            const float ac_end = 15.f;

            compression_curve.resize(n_moments);
            compression_curve[0] = compression_dc;

            for (size_t i = 1; i < n_moments; i++) {
                compression_curve[i] = deterministic_compression(compression_ac1, ac_end, n_moments, i);
            }

            timing = 0;
            break;

            break;
    }

    assert(compression_curve.size() == n_moments);
}


double linear_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort
    )
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    linear_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = linear_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = linear_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return linear_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double linavg_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort
    )
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    linavg_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = linavg_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = linavg_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return linavg_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double unbounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    unbounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = unbounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = unbounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return unbounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    bounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = bounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = bounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return bounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double unbounded_to_bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    unbounded_to_bounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = unbounded_to_bounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = unbounded_to_bounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return unbounded_to_bounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double upperbound_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_max;

    upperbound_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments,
        relative_scales,
        global_max
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = upperbound_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs,
        relative_scales,
        global_max
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = upperbound_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs,
                relative_scales,
                global_max
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return upperbound_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        relative_scales,
        global_max,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


double twobounds_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_min, global_max;

    twobounds_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        normalize_moments,
        relative_scales,
        global_min,
        global_max
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);

    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> q_normalized_moments;

    quantize_dequantize_image(
        normalized_moments,
        q_normalized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        q_normalized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        subsampling_factor_curve[1],
        1, compression_curve[1],
        effort
    );

    assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

    const double base_err = twobounds_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max
    );

    const float frame_distance_stop_inc = 0.1f;
    const float max_frame_distance = 15.0f;

#ifdef DEBLOG
    std::cout << "Starting compression curve computation..." << std::endl;
    std::cout << "start err: " << base_err << std::endl;
#endif // DEBLOG

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

#ifdef DEBLOG
        std::cout << "    INIT ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG

        float low = compression_curve[m];
        float high = max_frame_distance;
        float frame_distance = (low + high) / 2.f;

        while(low <= high) {
            compress_decompress_single_image(
                q_normalized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                subsampling_factor_curve[m],
                m, frame_distance,
                effort
            );

            assert(compressed_decompressed_moments.size() == q_normalized_moments.size());

            const double curr_err = twobounds_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs,
                relative_scales,
                global_min,
                global_max
            );

            float next_frame_distance;

            // Move left of right depending on the error
            if (base_err < curr_err) {
                // mv left
                high = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (next_frame_distance - low < frame_distance_stop_inc) {
                    frame_distance = low;
                    break;
                }
            } else {
                // mv right
                low = frame_distance;
                next_frame_distance = (low + high) / 2.f;

                if (high - next_frame_distance < frame_distance_stop_inc) {
                    break;
                }
            }

            frame_distance = next_frame_distance;

#ifdef DEBLOG
            std::cout << "           moment[" << m << "] d: " << frame_distance << " e: " << curr_err << std::endl;
#endif // DEBLOG
        }

        compression_curve[m] = frame_distance;

#ifdef DEBLOG
        std::cout << "    END  ~ moment[" << m << "] d: " << compression_curve[m] << std::endl;
#endif // DEBLOG
    }

    return twobounds_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );
}


/*****************************************************************************/
/* Error for a compression curve                                             */
/*****************************************************************************/

double linear_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    linear_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double linavg_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    linavg_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double unbounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    unbounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double bounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    bounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double unbounded_to_bounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    unbounded_to_bounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double upperbound_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    upperbound_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        relative_scales,
        global_max,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}


double twobounds_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_min,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);
    assert(subsampling_factor_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        subsampling_factor_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    twobounds_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return Util::rmse_images(
        ref_spectral_image,
        decompressed_spectral_image,
        width * height,
        wavelengths.size()
    );
}

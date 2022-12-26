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

#include "curve_compression.h"

#include <jxl/encode.h>
#include <jxl/encode_cxx.h>

#include <jxl/decode.h>
#include <jxl/decode_cxx.h>

#include <jxl/thread_parallel_runner.h>
#include <jxl/thread_parallel_runner_cxx.h>

#include "moments_image.h"
#include "moments_error.h"
#include "Util.h"

#include <vector>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <iostream>

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
    uint32_t bits_per_sample,
    uint32_t exponent_bits_per_sample,
    float frame_distance,
    uint32_t downsampling_ratio)
{
    assert(framebuffer_in.size() == width * height);
    framebuffer_out.resize(framebuffer_in.size());

    // 1. Compress the framebuffer using JXL
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

    JxlEncoderInitBasicInfo(&basic_info);

    basic_info.xsize                    = width;
    basic_info.ysize                    = height;
    basic_info.num_extra_channels       = 0;
    basic_info.num_color_channels       = 1;
    basic_info.bits_per_sample          = bits_per_sample;
    basic_info.exponent_bits_per_sample = exponent_bits_per_sample;
    basic_info.uses_original_profile    = JXL_TRUE;

    enc_status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);

    CHECK_JXL_ENC_STATUS(enc_status);

    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

    // Set compression quality
    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
    CHECK_JXL_ENC_STATUS(enc_status);

    if (frame_distance > 0) {
        enc_status = JxlEncoderSetFrameLossless(frame_settings, JXL_FALSE);
        enc_status = JxlEncoderSetFrameDistance(frame_settings, frame_distance);
    } else {
        enc_status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    }
    CHECK_JXL_ENC_STATUS(enc_status);

    enc_status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, downsampling_ratio);
    CHECK_JXL_ENC_STATUS(enc_status);

    color_encoding.color_space = JXL_COLOR_SPACE_GRAY;
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

    enc_status = JxlEncoderAddImageFrame(
        frame_settings,
        &format,
        framebuffer_in.data(),
        width * height * sizeof(float)
    );

    CHECK_JXL_ENC_STATUS(enc_status);

    JxlEncoderCloseInput(enc.get());

    std::vector<uint8_t> compressed(64);

    uint8_t* next_out = compressed.data();
    size_t avail_out = compressed.size() - (next_out - compressed.data());

    enc_status = JXL_ENC_NEED_MORE_OUTPUT;

    // for (status = JXL_ENC_NEED_MORE_OUTPUT ;; status == JXL_ENC_NEED_MORE_OUTPUT) {
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
    uint32_t bits_per_sample,
    size_t i, float frame_distance)
{
    assert(input_image.size() == width * height * n_moments);

    // Extract a single moment from the input image
    std::vector<float> input_framebuffer(width * height);
    std::vector<float> compressed_framebuffer;

    for (size_t px = 0; px < width * height; px++) {
        input_framebuffer[px] = (float)input_image[px * n_moments + i];
    }

    compress_decompress_framebuffer(
        input_framebuffer,
        compressed_framebuffer,
        width, height, bits_per_sample,
        0, frame_distance, 1
    );

    // Copy back the data to the output image
    output_image.resize(width * height * n_moments);

    std::memcpy(
        output_image.data(),
        input_image.data(),
        width * height * n_moments * sizeof(double)
    );

    // Change what has changed
    for (size_t px = 0; px < width * height; px++) {
        output_image[px * n_moments + i] = (double)compressed_framebuffer[px];
    }
}


void compress_decompress_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(input_image.size() == width * height * n_moments);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    output_image.resize(width * height * n_moments);

    // JXL compression for every layers
    #pragma omp parallel for
    for (size_t m = 0; m < n_moments; m++) {
        std::vector<float> framebuffer_in(width * height);
        std::vector<float> framebuffer_out;

        // Copy cast to float
        for (size_t px = 0; px < width * height; px++) {
            framebuffer_in[px] = input_image[px * n_moments + m];
        }

        // Compress / decompress
        const int bps = quantization_curve[m];
        int exponent_bits = 0;

        if (m == 0) {
            if (bps == 16) {
                exponent_bits = 5;
            } else if (bps == 32) {
                exponent_bits = 8;
            } else {
                throw std::runtime_error("Unknown quantization ratio for 0th moment");
            }
        }

        // TODO: double check that the quantization is applied
        compress_decompress_framebuffer(
            framebuffer_in,
            framebuffer_out,
            width, height,
            bps, exponent_bits,
            compression_curve[m],
            1);

        // Copy cast to float
        for (size_t px = 0; px < width * height; px++) {
            output_image[px * n_moments + m] = framebuffer_out[px];
        }
    }
}


/*****************************************************************************/
/* Create compression curves                                                 */
/*****************************************************************************/

double linear_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    linear_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

    const double base_err = linear_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

            const double curr_err = linear_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
    }

    return linear_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        compression_curve
    );
}


double unbounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    unbounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

    const double base_err = unbounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

            const double curr_err = unbounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
    }

    return unbounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        compression_curve
    );
}


double bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    bounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

    const double base_err = bounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

            const double curr_err = bounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
    }

    return bounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        compression_curve
    );
}


double unbounded_to_bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    unbounded_to_bounded_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

    const double base_err = unbounded_to_bounded_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

            const double curr_err = unbounded_to_bounded_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
    }

    return unbounded_to_bounded_error_for_compression_curve(
        wavelengths, spectral_image,
        width, height,
        n_moments,
        normalized_moments,
        mins, maxs,
        quantization_curve,
        compression_curve
    );
}


double upperbound_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_max;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    upperbound_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        relative_scales,
        global_max
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

    const double base_err = upperbound_average_err(
        wavelengths,
        spectral_image,
        width * height, n_moments,
        compressed_decompressed_moments,
        mins, maxs,
        relative_scales,
        global_max
    );

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

            const double curr_err = upperbound_average_err(
                wavelengths,
                spectral_image,
                width * height, n_moments,
                compressed_decompressed_moments,
                mins, maxs,
                relative_scales,
                global_max
            );

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
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
        compression_curve
    );
}


double twobounds_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve)
{
    assert(spectral_image.size() == width * height * wavelengths.size());
    assert(quantization_curve.size() == n_moments);

    compression_curve.resize(n_moments);

    std::vector<double> normalized_moments;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_min, global_max;

    const float frame_distance_inc = 0.1f;
    const float max_frame_distance = 15.0f;

    twobounds_compress_spectral_image(
        wavelengths, spectral_image,
        width * height, n_moments,
        normalized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max
    );

    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);

    // TODO: Check if this is really necessary, it is skipped in other parts of the code
    // Quantize / dequantize each moment based on the provided quantization
    // curve
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < width * height; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                Util::quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    compression_curve[0] = compression_dc;
    compression_curve[1] = compression_ac1;

    std::vector<double> compressed_decompressed_moments;

    compress_decompress_single_image(
        quantized_moments,
        compressed_decompressed_moments,
        width, height, n_moments,
        quantization_curve[1],
        1, compression_curve[1]
    );

    assert(compressed_decompressed_moments.size() == quantized_moments.size());

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

    for (size_t m = 2; m < n_moments; m++) {
        compression_curve[m] = compression_curve[m - 1];

        for (float frame_distance = compression_curve[m] + frame_distance_inc; frame_distance <= max_frame_distance; frame_distance += frame_distance_inc) {
            compress_decompress_single_image(
                quantized_moments,
                compressed_decompressed_moments,
                width, height, n_moments,
                quantization_curve[m],
                m, frame_distance
            );

            assert(compressed_decompressed_moments.size() == quantized_moments.size());

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

            if (curr_err >= base_err) {
                break;
            }

            compression_curve[m] = frame_distance;
        }
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
        compression_curve
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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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


double unbounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
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

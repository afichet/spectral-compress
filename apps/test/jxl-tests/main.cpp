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
#include <vector>
#include <sstream>
#include <cmath>

#include <JXLImage.h>
#include <Util.h>
#include <curve_quantization.h>

#include <curve_compression.h>

#include <sstream>


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

void compress_framebuffer(
    const std::vector<float>& framebuffer_in,
    const char* filename,
    uint32_t width, uint32_t height,
    uint32_t bits_per_sample,
    uint32_t exponent_bits_per_sample,
    float frame_distance,
    uint32_t subsampling_ratio)
{
    assert(framebuffer_in.size() == width * height);

    // 1. Compress the framebuffer using JXL
    JxlEncoderStatus           status;
    JxlEncoderPtr              enc;
    JxlThreadParallelRunnerPtr runner;
    JxlBasicInfo               basic_info;
    JxlColorEncoding           color_encoding;

    runner = JxlThreadParallelRunnerMake(
            nullptr,
            JxlThreadParallelRunnerDefaultNumWorkerThreads()
        );

    enc = JxlEncoderMake(nullptr);

    status = JxlEncoderSetParallelRunner(
        enc.get(),
        JxlThreadParallelRunner,
        runner.get()
    );

    CHECK_JXL_ENC_STATUS(status);

    JxlEncoderInitBasicInfo(&basic_info);

    basic_info.xsize                    = width;
    basic_info.ysize                    = height;
    basic_info.num_extra_channels       = 0;
    basic_info.num_color_channels       = 1;
    basic_info.bits_per_sample          = bits_per_sample;
    basic_info.exponent_bits_per_sample = exponent_bits_per_sample;
    basic_info.uses_original_profile    = JXL_TRUE;

    status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------

    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

    // Set compression quality
    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
    CHECK_JXL_ENC_STATUS(status);

    if (frame_distance > 0) {
        // FIXME: quantization is ignored
        status = JxlEncoderSetFrameLossless(frame_settings, JXL_FALSE);
        status = JxlEncoderSetFrameDistance(frame_settings, frame_distance);
    } else {
        status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    }
    CHECK_JXL_ENC_STATUS(status);

    status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, subsampling_ratio);
    CHECK_JXL_ENC_STATUS(status);


    color_encoding.color_space       = JXL_COLOR_SPACE_GRAY;
    color_encoding.white_point       = JXL_WHITE_POINT_D65;
    color_encoding.primaries         = JXL_PRIMARIES_SRGB;
    color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
    color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

    status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);

    CHECK_JXL_ENC_STATUS(status);

    // JxlPixelFormat format;
    // format.num_channels = 1;
    // format.data_type    = JXL_TYPE_FLOAT,
    // format.endianness   = JXL_NATIVE_ENDIAN;
    // format.align        = 0;

    // status = JxlEncoderAddImageFrame(
    //     frame_settings,
    //     &format,
    //     framebuffer_in.data(),
    //     width * height * sizeof(float)
    // );

    std::vector<uint8_t> converted_framebuffer(width * height);
    for (size_t i = 0; i < width * height; i++) {
        converted_framebuffer[i] = 255.f * framebuffer_in[i];
    }
    JxlPixelFormat format;
    format.num_channels = 1;
    format.data_type    = JXL_TYPE_UINT8,
    format.endianness   = JXL_NATIVE_ENDIAN;
    format.align        = 0;

    status = JxlEncoderAddImageFrame(
        frame_settings,
        &format,
        converted_framebuffer.data(),
        width * height * sizeof(uint8_t)
    );


    CHECK_JXL_ENC_STATUS(status);

    JxlEncoderCloseInput(enc.get());

    std::vector<uint8_t> compressed(64);

    uint8_t* next_out = compressed.data();
    size_t avail_out = compressed.size() - (next_out - compressed.data());

    status = JXL_ENC_NEED_MORE_OUTPUT;

    while (status == JXL_ENC_NEED_MORE_OUTPUT) {
        status = JxlEncoderProcessOutput(enc.get(), &next_out, &avail_out);

        if (status == JXL_ENC_NEED_MORE_OUTPUT) {
            size_t offset = next_out - compressed.data();
            compressed.resize(compressed.size() * 2);

            next_out  = compressed.data() + offset;
            avail_out = compressed.size() - offset;
        }
    }

    CHECK_JXL_ENC_STATUS(status);

    compressed.resize(next_out - compressed.data());
    // TODO: Remove
    // Write file
    std::FILE* file = std::fopen(filename, "wb");
    std::fwrite(compressed.data(), sizeof(uint8_t), compressed.size(), file);
    std::fclose(file);
}

/**
 * Performs various tests on the JXL compressor
 * - Quantization: we expect progressively high banding
 * - Subsampling: we expect increasingly lower resolution
 */
int main(int argc, char* argv[])
{
    (void)argc; (void)argv;

    size_t width = 640;
    size_t height = 480;

    std::vector<float> framebuffer(width * height);

    // ------------------------------------------------------------------------
    // Generate a gradiant to check if the quantization ratio works as expected
    // ------------------------------------------------------------------------

    for (size_t y = 0; y < height; y++) {
        const float value = (float)y / (float)(height - 1);

        for (size_t x = 0; x < width; x++) {
            framebuffer[y * width + x] = value;
        }
    }

    // Save multiple files with a different quantization
    for (size_t q = 8; q > 1; q--) {
        std::stringstream ss;
        ss << "quantization_" << q << ".jxl";

        JXLImage jxl_out(width, height);
        jxl_out.appendFramebuffer(
            framebuffer, 1,
            std::make_pair(q, 0),
            1, 0.1f
        );

        jxl_out.write(ss.str());
    }

    // 2nd attempt in a single pass
    JXLImage jxl_out(width, height);

    for (size_t q = 8; q > 1; q--) {
        jxl_out.appendFramebuffer(
            framebuffer, 1,
            std::make_pair(q, 0),
            1, 0.1f
        );
    }

    jxl_out.write("jxl_quantization.jxl");

    // With compress decompress
    for (size_t q = 8; q > 1; q--) {
        std::stringstream ss;
        ss << "quantization_cc_8_" << q << ".jxl";

        compress_framebuffer(
            framebuffer,
            ss.str().c_str(),
            width, height,
            q, 0,
            0, 1
        );
    }
    // std::vector<float> foo;
    // compress_decompress_framebuffer(
    //     framebuffer, foo,
    //     width, height, 3, 0, 0.1, 1);


    // Save multiple files with a different quantization using the
    // "sotfware" quantization function
    std::vector<double> framebuffer_d;
    Util::cast_vector(framebuffer, framebuffer_d);

    for (size_t q = 8; q > 1; q--) {
        std::vector<double> qdq_d;
        std::vector<float> qdq;

        quantize_dequantize_single_image(
            framebuffer_d, qdq_d,
            width * height,
            1,
            0, q
        );

        Util::cast_vector(qdq_d, qdq);

        std::stringstream ss;
        ss << "sw_quantization_" << q << ".jxl";

        JXLImage jxl_out(width, height);
        jxl_out.appendFramebuffer(
            qdq, 1,
            std::make_pair(32, 8),
            1, 0.f
        );

        jxl_out.write(ss.str());
    }

    // ------------------------------------------------------------------------
    // Generate lines of different heigth to check if the subsampling works
    // as expected
    // ------------------------------------------------------------------------

    int start_sz = 1;
    int remaining_lines = start_sz;
    int curr_sz = start_sz;

    for (size_t y = 0; y < height; y++) {
        float value = 0;

        if (remaining_lines > 0) {
            value = 1;
        } else if (-remaining_lines < curr_sz) {
            value = 0;
        } else {
            curr_sz += 1;
            remaining_lines = curr_sz;
        }

        for (size_t x = 0; x < width; x++) {
            framebuffer[y * width + x] = value;
        }

        --remaining_lines;
    }

    // It seems JXL supports this subsampling rates only
    size_t subsampling_ratios[] = {1, 2, 4, 8};

    for (size_t subsampling: subsampling_ratios) {
        std::stringstream ss;
        ss << "downsampled_" << subsampling << ".jxl";

        JXLImage jxl_out(width, height);
        jxl_out.appendFramebuffer(
            framebuffer, 1,
            std::make_pair(8, 0),
            subsampling, 0.f
        );

        jxl_out.write(ss.str());
    }

    return 0;
}

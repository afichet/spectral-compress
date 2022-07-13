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

#pragma once

#include <cstdint>
#include <vector>
#include <string>

#include <jxl/encode.h>
#include <jxl/encode_cxx.h>

#include <jxl/decode.h>
#include <jxl/decode_cxx.h>

#include <jxl/thread_parallel_runner.h>
#include <jxl/thread_parallel_runner_cxx.h>

#include <SGEG_box.h>

class JXLImageWriter
{
public:
    JXLImageWriter(
        size_t width,
        size_t height,
        SGEG_box sgeg_box,
        size_t n_sub_framebuffers = 0);

    void addMainFramebuffer(
        void *framebuffer,
        JxlPixelFormat pixel_format,
        size_t el_size,
        int32_t downsampling = 1);

    void addMainFramebuffer(float *framebuffer, int32_t downsampling = 1);
    void addMainFramebuffer(uint8_t *framebuffer, int32_t downsampling = 1);
    void addMainFramebuffer(uint16_t *framebuffer, int32_t downsampling = 1);

    void addSubFramebuffer(
        void *framebuffer,
        JxlPixelFormat pixel_format,
        JxlExtraChannelInfo extra_info,
        size_t el_size,
        size_t index,
        int32_t downsampling = 1,
        const char *channel_name = nullptr);

    void addSubFramebuffer(float *framebuffer, size_t index, int32_t downsampling = 1, const char *channel_name = nullptr);
    void addSubFramebuffer(uint8_t *framebuffer, size_t index, int32_t downsampling = 1, const char *channel_name = nullptr);
    void addSubFramebuffer(uint16_t *framebuffer, size_t index, int32_t downsampling = 1, const char *channel_name = nullptr);

    void save(const char *filename);

private:
    JxlEncoderPtr              _enc;
    JxlThreadParallelRunnerPtr _runner;
    JxlBasicInfo               _basic_info;
    JxlColorEncoding           _color_encoding;
    SGEG_box                   _sgeg_box;
};


class JXLImageReader
{
public:
    JXLImageReader(const char *filename);
    virtual ~JXLImageReader();

    uint32_t width() const;
    uint32_t height() const;
    uint32_t n_subframebuffers() const;

    SGEG_box get_sgeg() const;

    void getMainFramebuffer(std::vector<float> &framebuffer) const;
    void getSubFramebuffer(std::vector<float> &framebuffer, size_t index) const;

    void print_basic_info() const;

private:
    JxlDecoderPtr              _dec;
    JxlThreadParallelRunnerPtr _runner;
    JxlBasicInfo               _basic_info;
    SGEG_box                   _sgeg_box;

    std::vector<float>              _main_framebuffer;
    std::vector<std::vector<float>> _sub_framebuffers;
    std::vector<std::string>        _sub_framebuffers_names;
};
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


class JXLFramebuffer
{
public:
    JXLFramebuffer(
        uint32_t width, uint32_t height,
        uint32_t n_color_channels = 1,
        uint32_t n_bits_per_sample = 32, 
        uint32_t n_exponent_bits_per_sample = 8,
        uint32_t downsampling_factor = 1,
        const char* name = nullptr);

    JXLFramebuffer(
        const std::vector<float>& framebuffer,
        uint32_t n_color_channels = 1,
        uint32_t n_bits_per_sample = 32, 
        uint32_t n_exponent_bits_per_sample = 8,
        uint32_t downsampling_factor = 1,
        const char* name = nullptr);

    virtual ~JXLFramebuffer();

    uint32_t getNColorChannels() const { return _pixel_format.num_channels; }
    uint32_t getBitsPerSample()  const { return _n_bits_per_sample; }
    uint32_t getExponentBitsPerSample() const { return _n_exponent_bits_per_sample; }
    uint32_t getDownsamplingFactor() const { return _downsampling_factor; }

    JxlPixelFormat getPixelFormat() const { return _pixel_format; }

    const char* getName() const { return _name; }

    void setName(const char* name);

    size_t getSizeBytes() const {
        return _pixel_data.size() * sizeof(float);
    }

protected:
    char* _name;

    uint32_t _n_bits_per_sample;
    uint32_t _n_exponent_bits_per_sample;
    uint32_t _downsampling_factor; // Not possible with the current state of the API

    JxlPixelFormat _pixel_format;

public:
    std::vector<float> _pixel_data;
};



class JXLImage
{
public:
    JXLImage(uint32_t width, uint32_t height);
    JXLImage(const char* filename);

    virtual ~JXLImage();

    void appendFramebuffer(
        const std::vector<float>& framebuffer,
        uint32_t n_channels,
        uint32_t enc_bits_per_sample = 32,
        uint32_t enc_exponent_bits_per_sample = 8,
        uint32_t enc_downsampling_factor = 1,
        const char* name = nullptr);

    void write(const char* filename) const;

    void setBox(const SGEG_box& box);
    SGEG_box getBox() const { return _sgeg_box; }

    uint32_t width()  const { return _width; }
    uint32_t height() const { return _height; }

    size_t n_framebuffers() const { return _framebuffers.size(); }
    
    std::vector<float>& getFramebufferData(size_t index) const
    {
        return _framebuffers[index]->_pixel_data;
    }

    JXLFramebuffer* getFramebuffer(size_t index) const
    {
        return _framebuffers[index];
    }

protected:
    uint32_t _width;
    uint32_t _height;

    SGEG_box _sgeg_box;

    std::vector<JXLFramebuffer*> _framebuffers;
};
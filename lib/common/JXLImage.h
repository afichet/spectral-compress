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

#include <SGEGBox.h>


class JXLFramebuffer
{
public:
    JXLFramebuffer(
        uint32_t width, uint32_t height,
        uint32_t n_color_channels = 1,
        uint32_t n_bits_per_sample = 32,
        uint32_t n_exponent_bits_per_sample = 8,
        uint32_t downsampling_factor = 1,
        float    framedistance = 0,
        const char* name = nullptr);

    JXLFramebuffer(
        const std::vector<float>& framebuffer,
        uint32_t n_color_channels = 1,
        uint32_t n_bits_per_sample = 32,
        uint32_t n_exponent_bits_per_sample = 8,
        uint32_t downsampling_factor = 1,
        float    framedistance = 0,
        const char* name = nullptr);

    virtual ~JXLFramebuffer();

    uint32_t getNColorChannels() const { return _pixel_format.num_channels; }
    uint32_t getBitsPerSample()  const { return _n_bits_per_sample; }
    uint32_t getExponentBitsPerSample() const { return _n_exponent_bits_per_sample; }
    uint32_t getDownsamplingFactor() const { return _downsampling_factor; }
    float getFramedistance() const { return _framedistance; }

    JxlPixelFormat getPixelFormat() const { return _pixel_format; }

    std::vector<float>& getPixelData() { return _pixel_data; }
    const std::vector<float>& getPixelDataConst() const { return _pixel_data; }

    const char* getName() const { return _name; }

    void setName(const char* name);

    size_t getSizeBytes() const {
        return _pixel_data.size() * sizeof(float);
    }

    void dump(FILE* stream) const;

    static JXLFramebuffer* read_dump(FILE* stream);

protected:
    char* _name;

    uint32_t _n_bits_per_sample;
    uint32_t _n_exponent_bits_per_sample;
    uint32_t _downsampling_factor; // Not possible with the current state of the API
    float    _framedistance;

    JxlPixelFormat _pixel_format;

    std::vector<float> _pixel_data;
};


class JXLImage
{
public:
    JXLImage(const char* filename);
    JXLImage(const std::string& filename);

    JXLImage(uint32_t width, uint32_t height);

    virtual ~JXLImage();

    size_t appendFramebuffer(
        const std::vector<float>& framebuffer,
        uint32_t n_channels,
        uint32_t enc_bits_per_sample = 32,
        uint32_t enc_exponent_bits_per_sample = 8,
        uint32_t enc_downsampling_factor = 1,
        float    enc_framedistance = .1f,
        const char* name = nullptr);

    void write(const char* filename) const;
    void write(const std::string& filename) const;

    void setBox(const SGEGBox& box);
    const SGEGBox& getBox() const { return _sgeg_box; }

    uint32_t width()  const { return _width; }
    uint32_t height() const { return _height; }

    size_t n_framebuffers() const { return _framebuffers.size(); }

    const std::vector<float>& getFramebufferDataConst(size_t index) const
    {
        return _framebuffers[index]->getPixelDataConst();
    }

    std::vector<float>& getFramebufferData(size_t index) const
    {
        return _framebuffers[index]->getPixelData();
    }

    JXLFramebuffer* getFramebuffer(size_t index) const
    {
        return _framebuffers[index];
    }

    std::vector<JXLFramebuffer*>& getFramebuffers() {
        return _framebuffers;
    }

    const std::vector<JXLFramebuffer*>& getFramebuffersConst() const {
        return _framebuffers;
    }

    void dump(const char* filename) const;
    void dump(const std::string& filename) const;

    static JXLImage* read_dump(const char* filename);

protected:
    void load(const char* filename);

    uint32_t _width;
    uint32_t _height;
    uint32_t _n_parts;

    SGEGBox _sgeg_box;

    std::vector<JXLFramebuffer*> _framebuffers;
};

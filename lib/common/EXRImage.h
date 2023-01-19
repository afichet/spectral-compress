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

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

#include <OpenEXR/ImfIO.h>
#include <OpenEXR/ImfCompression.h>


class EXRFramebuffer
{
public:
    EXRFramebuffer(
        uint32_t width, uint32_t height,
        const char* name);

    EXRFramebuffer(
        const std::vector<float>& framebuffer,
        const char* name);

    virtual ~EXRFramebuffer();

    const char* getName() const { return _name.data(); }
    const std::vector<float>& getPixelDataConst() const { return _pixel_data; }
    std::vector<float>& getPixelData() { return _pixel_data; }

protected:
    std::vector<float> _pixel_data;
    std::vector<char> _name;
};


class EXRImage
{
public:
    EXRImage(const char *filename);

    EXRImage(const std::string& filename);

    EXRImage(uint32_t width, uint32_t height, Imf::Compression compression = Imf::Compression::ZIP_COMPRESSION);

    virtual ~EXRImage();

    void appendFramebuffer(
        const std::vector<float>& framenuffer,
        const char* name);

    void write(const char* filename) const;

    void write(const std::string& filename) const;

    uint32_t width()  const { return _width; }
    uint32_t height() const { return _height; }

    size_t n_framebuffers() const { return _framebuffers.size(); }

    std::vector<EXRFramebuffer*>& getFramebuffers() {
        return _framebuffers;
    }

    const std::vector<EXRFramebuffer*>& getFramebuffersConst() const {
        return _framebuffers;
    }

    void setAttributesData(const std::vector<uint8_t>& data) { _attributes_data = data; }

    const std::vector<uint8_t>& getAttributesData() const { return _attributes_data; }

protected:
    void load(const char* filename);

    uint32_t _width, _height;
    Imf::Compression _compression;
    std::vector<EXRFramebuffer*> _framebuffers;

    std::vector<uint8_t> _attributes_data;
};

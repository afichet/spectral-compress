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
#include <cstddef>
#include <vector>
#include <string>


enum PixelType {
    UINT = 0,
    HALF = 1,
    FLOAT = 2,
    N_PIXEL_TYPES
};

struct SpectralFramebuffer
{
    std::string root_name;
    std::vector<float> wavelengths_nm;
    std::vector<float> image_data;

    PixelType pixel_type;
};


struct GreyFramebuffer
{
    std::string layer_name;
    std::vector<float> image_data;

    PixelType pixel_type;
};


class EXRSpectralImage
{
public:
    EXRSpectralImage(const char *filename);

    EXRSpectralImage(uint32_t width, uint32_t height);

    virtual ~EXRSpectralImage();

    void appendSpectralFramebuffer(
        const std::vector<float> &wavelengths_nm,
        const std::vector<float> &framebuffer,
        const std::string& prefix,
        const PixelType save_as_type = PixelType::FLOAT);

    void appendSpectralFramebuffer(
        const std::vector<float> &wavelengths_nm,
        const std::vector<float> &framebuffer,
        const char* prefix,
        const PixelType save_as_type = PixelType::FLOAT);

    void appendExtraFramebuffer(
        const std::vector<float>& framenuffer,
        const std::string& name,
        const PixelType save_as_type = PixelType::FLOAT);

    void appendExtraFramebuffer(
        const std::vector<float>& framenuffer,
        const char* name,
        const PixelType save_as_type = PixelType::FLOAT);

    void write(const char* filename) const;

    uint32_t width()  const;
    uint32_t height() const;

    std::vector<SpectralFramebuffer*>& getSpectralFramebuffers();

    std::vector<GreyFramebuffer*>& getExtraFramebuffers();

    void setAttributesData(const std::vector<uint8_t>& data);

    const std::vector<uint8_t>& getAttributesData() const;

    static double to_nm(
        const std::string& value,
        const std::string& prefix,
        const std::string& units);

    /**
     * @brief Dumps the data to a file
     *
     * This produces a raw file, not an OpenEXR, with a dump of the content of
     * the class
     *
     * @param filename Filename to dump the data to
     */
    void dump(const char* filename) const;

    static EXRSpectralImage* read_dump(const char* filename);

protected:
    uint32_t _width, _height;

    std::vector<SpectralFramebuffer*> _spectral_framebuffers;
    std::vector<GreyFramebuffer*> _extra_framebuffers;

    std::vector<uint8_t> _attributes_data;
};
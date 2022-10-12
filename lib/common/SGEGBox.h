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

#include <vector>
#include <cstddef>
#include <cstdint>
#include <string>


struct SGEGSpectralGroup
{
    std::vector<char> root_name;
    std::vector<uint32_t> layer_indices;
    std::vector<float> wavelengths;
    std::vector<float> mins;
    std::vector<float> maxs;    

    size_t getRaw(uint8_t data[]) const;
    size_t fromRaw(const uint8_t data[]);

    size_t size() const;
    SGEGSpectralGroup& operator=(const SGEGSpectralGroup& other);
};


struct SGEGGrayGroup
{
    std::vector<char> layer_name;
    uint32_t layer_index;

    size_t getRaw(uint8_t data[]) const;
    size_t fromRaw(const uint8_t data[]);

    size_t size() const;
    SGEGGrayGroup& operator=(const SGEGGrayGroup& other);
};


// Spectral Graphics Extended Group

struct SGEGBox
{
    uint32_t revision;
    uint32_t n_parts;

    std::vector<SGEGSpectralGroup> spectral_groups;
    std::vector<SGEGGrayGroup> gray_groups;

    std::vector<uint8_t> exr_attributes;

    SGEGBox();

    SGEGBox(const std::vector<uint8_t> &data);

    SGEGBox(const SGEGBox& other);

    virtual ~SGEGBox();

    void getRaw(std::vector<uint8_t> &data) const;

    size_t getRaw(uint8_t data[]) const;

    size_t size() const;

    SGEGBox& operator=(const SGEGBox& other);

    void print() const;
};

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

// class ImageMeta
// {
// public:
//     uint32_t n_moments;
//     SpectrumType spectrum_type;
//     uint32_t len_channel_prefix;
//     char* channel_prefix;
// }


// Spectral Gamma Extended Graphics

class SGEG_box
{
public:
    uint32_t n_moments;
    uint32_t n_wavelengths;

    std::vector<float> moment_min;
    std::vector<float> moment_max;
    std::vector<float> wavelengths;

    std::vector<uint32_t> subimage_idx;

    bool is_reflective;

    SGEG_box(uint32_t n_moments = 0, uint32_t n_wavelengths = 0);

    SGEG_box(const std::vector<uint8_t> &data);

    SGEG_box(const SGEG_box& other);

    void getRaw(std::vector<uint8_t> &data) const;

    void print_info() const;

    SGEG_box& operator=(const SGEG_box& other)
    {
        if (this != &other) {
            n_moments     = other.n_moments;
            n_wavelengths = other.n_wavelengths;
            
            moment_min    = other.moment_min;
            moment_max    = other.moment_max;
            wavelengths   = other.wavelengths;
            
            subimage_idx  = other.subimage_idx;

            is_reflective = other.is_reflective;
        }

        return *this;
    }

    size_t size() const;
};

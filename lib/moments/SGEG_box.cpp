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

#include "SGEG_box.h"

#include <cstring>
#include <cstdint>


SGEG_box::SGEG_box(uint32_t sz)
{
    n_moments = sz;
    moment_min.resize(n_moments);
    moment_max.resize(n_moments);
}


SGEG_box::SGEG_box(const std::vector<uint8_t> &data)
{
    size_t offset = 0;
    size_t next_offset = 0;

    // ------------------------------------------------------------------------
    // Number of moments
    // ------------------------------------------------------------------------
    next_offset = sizeof(uint32_t);

    memcpy(&n_moments, data.data(), next_offset);
    offset += next_offset;

    moment_min.resize(n_moments);
    moment_max.resize(n_moments);

    // ------------------------------------------------------------------------
    // Rescaling factor [min..max]
    // ------------------------------------------------------------------------

    next_offset = n_moments * sizeof(float);

    memcpy(moment_min.data(), data.data() + offset, next_offset);
    offset += next_offset;

    memcpy(moment_max.data(), data.data() + offset, next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Original spectral bounds [\lambda_min..\lambda_max]
    // ------------------------------------------------------------------------

    next_offset = sizeof(float);

    memcpy(&wl_min_nm, data.data() + offset, next_offset);
    offset += next_offset;

    memcpy(&wl_max_nm, data.data() + offset, next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Original number of spectral samples
    // ------------------------------------------------------------------------

    next_offset = sizeof(uint32_t);

    memcpy(&n_wl_original, data.data() + offset, next_offset);
}


SGEG_box::SGEG_box(const SGEG_box& other)
    : n_moments(other.n_moments)
    , moment_min(other.moment_min)
    , moment_max(other.moment_max)
    , wl_min_nm(other.wl_min_nm)
    , wl_max_nm(other.wl_max_nm)
    , n_wl_original(other.n_wl_original)
{
}


void SGEG_box::getRaw(std::vector<uint8_t> &data) const
{
    size_t offset = 0;
    size_t next_offset = 0;

    data.resize(size());

    // ------------------------------------------------------------------------
    // Number of moments
    // ------------------------------------------------------------------------

    next_offset = sizeof(uint32_t);

    memcpy(data.data() + offset, &n_moments, next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Rescaling factor [min..max]
    // ------------------------------------------------------------------------

    next_offset = n_moments * sizeof(float);

    memcpy(data.data() + offset, moment_min.data(), next_offset);
    offset += next_offset;

    memcpy(data.data() + offset, moment_max.data(), next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Original spectral bounds [\lambda_min..\lambda_max]
    // ------------------------------------------------------------------------

    next_offset = sizeof(float);

    memcpy(data.data() + offset, &wl_min_nm, next_offset);
    offset += next_offset;

    memcpy(data.data() + offset, &wl_max_nm, next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Original number of spectral samples
    // ------------------------------------------------------------------------

    next_offset = sizeof(uint32_t);

    memcpy(data.data() + offset, &n_wl_original, next_offset);
}


size_t SGEG_box::size() const
{
    return sizeof(uint32_t) +              // n_moments
            n_moments * 2 * sizeof(float) + // [min..max]
            2 * sizeof(float) +             // [\lambda_min..\lambda_max]
            sizeof(uint32_t);               // n_spectral_samples
}
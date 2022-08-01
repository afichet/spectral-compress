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


SGEG_box::SGEG_box(uint32_t _n_moments, uint32_t _n_wavelengths)
    : n_moments(_n_moments)
    , n_wavelengths(_n_wavelengths)
    , moment_min(_n_moments)
    , moment_max(_n_moments)
    , wavelengths(_n_wavelengths)
{
}


SGEG_box::SGEG_box(const std::vector<uint8_t> &data)
{
    size_t offset = 0;
    size_t next_offset = 0;

    // ------------------------------------------------------------------------
    // Number of moments
    // ------------------------------------------------------------------------
    
    next_offset = sizeof(uint32_t);

    memcpy(&n_moments, data.data() + offset, next_offset);
    offset += next_offset;

    moment_min.resize(n_moments);
    moment_max.resize(n_moments);

    // ------------------------------------------------------------------------
    // Original number of wavelengths
    // ------------------------------------------------------------------------

    next_offset = sizeof(uint32_t);

    memcpy(&n_wavelengths, data.data() + offset, next_offset);
    offset += next_offset;

    wavelengths.resize(n_wavelengths);

    // ------------------------------------------------------------------------
    // Rescaling factor [min..max]
    // ------------------------------------------------------------------------

    next_offset = n_moments * sizeof(float);

    memcpy(moment_min.data(), data.data() + offset, next_offset);
    offset += next_offset;

    memcpy(moment_max.data(), data.data() + offset, next_offset);
    offset += next_offset;

    // ------------------------------------------------------------------------
    // Original spectral vlues
    // ------------------------------------------------------------------------

    next_offset = n_wavelengths * sizeof(float);

    memcpy(wavelengths.data(), data.data() + offset, next_offset);
    offset += next_offset;
}


SGEG_box::SGEG_box(const SGEG_box& other)
    : n_moments(other.n_moments)
    , n_wavelengths(other.n_wavelengths)
    , moment_min(other.moment_min)
    , moment_max(other.moment_max)
    , wavelengths(other.wavelengths)
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
    // Original number of wavelengths
    // ------------------------------------------------------------------------

    next_offset = sizeof(uint32_t);

    memcpy(data.data() + offset, &n_wavelengths, next_offset);
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
    // Original spectral values
    // ------------------------------------------------------------------------

    next_offset = n_wavelengths * sizeof(float);

    memcpy(data.data() + offset, wavelengths.data(), next_offset);
    offset += next_offset;
}


size_t SGEG_box::size() const
{
    return
        sizeof(uint32_t) +              // n_moments
        sizeof(uint32_t) +              // n_wavelengths
        2 * n_moments * sizeof(float) + // [min..max]
        n_wavelengths * sizeof(float);  // [\lambda_min..\lambda_max]
}


#include <iostream>

void SGEG_box::print_info() const
{
    std::cout << "    # moments: " << n_moments << std::endl;
    std::cout << "# wavelengths: " << n_wavelengths << std::endl;
    
    for (uint32_t i = 0; i < n_moments; i++) {
        std::cout << " + moment[" << i << "]"
                  <<  " - min, max: [" << moment_min[i] << ", " << moment_max[i] << "]" << std::endl;
    }
    
    for (uint32_t i = 0; i < n_wavelengths; i++) {
        std::cout << " + wavelength[" << i << "]"
                  << " - " << wavelengths[i] << " nm" << std::endl;
    }
}
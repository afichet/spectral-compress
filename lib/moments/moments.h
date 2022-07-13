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

#include <cstddef>
#include <vector>

// C interface

void wavelengths_to_phases(
    const float wavelengths[], 
    float phases[],
    size_t n_wavelengths);


void compute_moments(
    const float phases[],
    size_t n_phases,
    const float signal[],
    size_t n_moments, 
    float moments[]);


void solve_levinson(
    const float first_column[],
    float solution[],
    size_t size);


void dot_levinson(
    const float first_column[],
    float dot_product[],
    size_t size);


void levinson_from_dot(
    const float dot_product[],
    float first_column[],
    size_t size);


void compute_density(
    const float phases[],
    size_t n_phases,
    const float moments[],
    size_t n_moments,
    float density[]);


void compress_moments(
    const float moments[],
    size_t n_moments,
    float compressed_moments[]);


void decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[]);

// C++ Interface

void wavelengths_to_phases(
    const std::vector<float>& wavelengths, 
    std::vector<float>& phases);


void compute_moments(
    const std::vector<float>& phases, 
    const std::vector<float>& signal, 
    size_t n_moments,
    std::vector<float>& moments);


void solve_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& solution);


void dot_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& dot_product);


void levinson_from_dot(
    const std::vector<float>& dot_product,
    std::vector<float>& first_column);


void compute_density(
    const std::vector<float>& phases, 
    const std::vector<float>& moments,
    std::vector<float>& density);


void compress_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments);


void decompress_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments);

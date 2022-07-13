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

#include "moments.h"

#include <cmath>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <vector>
#include <complex>

// C interface

void wavelengths_to_phases(
    const float wavelengths[], 
    float phases[],
    size_t n_wavelengths)
{
    float min_wl = wavelengths[0];
    float max_wl = wavelengths[n_wavelengths - 1];

    for (size_t i = 0; i < n_wavelengths; i++) {
        phases[i] = M_PI * (wavelengths[i] - min_wl) / (max_wl - min_wl) - M_PI;
    }
}


void compute_moments(
    const float phases[],
    size_t n_phases,
    const float signal[],
    size_t n_moments, 
    float moments[])
{
    const float r = 1.f / (2.f * M_PI);

    for (size_t k = 0; k < n_moments + 1; k++) {
        moments[k] = 0;

        for (size_t p = 0; p < n_phases - 1; p++) {
            const float x_0 = phases[p + 0];
            const float x_1 = phases[p + 1];

            const std::complex<float> e_0 = std::polar(r, k * x_0);
            const std::complex<float> e_1 = std::polar(r, k * x_1);
            
            const float s_0 = (e_0 * signal[p + 0]).real();
            const float s_1 = (e_1 * signal[p + 1]).real();

            moments[k] += (x_1 - x_0) * (s_0 + s_1);
        }
    }
}


void solve_levinson(
    const float first_column[],
    float solution[],
    size_t size)
{
    solution[0] = 1.f / first_column[0];

    std::vector<float> temp_s(size);

    for (size_t i = 1; i < size; i++) {
        solution[i] = 0;

        float dot_product = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product * solution[i - k]) / (1.f - dot_product * dot_product);
        }

        memcpy(solution, temp_s.data(), temp_s.size() * sizeof(float));
    }
}


void dot_levinson(
    const float first_column[],
    float dot_product[],
    size_t size)
{
    std::vector<float> solution(size);
    std::vector<float> temp_s(size);

    dot_product[0] = 0;
    solution[0]    = 1.f / first_column[0];

    for (size_t i = 1; i < size; i++) {
        dot_product[i] = 0;
        solution[i]    = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product[i] += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * solution[i - k]) / (1.f - dot_product[i] * dot_product[i]);
        }

        solution = temp_s;
    }
}


void levinson_from_dot(
    const float dot_product[],
    float first_column[],
    size_t size)
{
    std::vector<float> solution(size);
    std::vector<float> temp_s(size);

    first_column[0] = dot_product[0];
    solution[0]     = 1.f / dot_product[0];

    for (size_t i = 1; i < size; i++) {
        const float radius = 1.f / solution[0]; // real??
        float center = 0;

        for (size_t k = 1; k < i; k++) {
            center += solution[k] * first_column[i - k];
        }

        first_column[i] = radius * (dot_product[i] - center);

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * solution[i - k]) / (1.f - dot_product[i] * dot_product[i]);
        }

        solution = temp_s;
    }
}


void compute_density(
    const float phases[],
    size_t n_phases,
    const float moments[],
    size_t n_moments,
    float density[])
{
    std::vector<float> v(n_moments + 1);
    
    for (size_t i = 0; i < v.size(); i++) {
        v[i] = moments[i] / (2.f * M_PI);
    }

    std::vector<float> solution;
    solve_levinson(v, solution);

    std::vector<std::complex<float>> denum(n_phases);

    const float r = 1.f / (2.f * M_PI);

    for (size_t p = 0; p < n_phases; p++) {
        std::complex<float> denum = 0;

        for (size_t k = 0; k < n_moments + 1; k++) {
            denum += solution[k] * std::polar(r, k * phases[p]);
        }

        denum = std::abs(denum);

        density[p] = r * solution[0] / (denum * denum).real();
    }
}


void compress_moments(
    const float moments[],
    size_t n_moments,
    float compressed_moments[])
{
    dot_levinson(moments, compressed_moments, n_moments + 1);
    compressed_moments[0] = moments[0];
}


void decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[])
{
    levinson_from_dot(compressed_moments, moments, n_compressed_moments + 1);
}


// C++ interface

void wavelengths_to_phases(
    const std::vector<float>& wavelengths, 
    std::vector<float>& phases) 
{
    phases.resize(wavelengths.size());
    wavelengths_to_phases(wavelengths.data(), phases.data(), wavelengths.size());
}


void compute_moments(
    const std::vector<float>& phases, 
    const std::vector<float>& signal, 
    size_t n_moments, 
    std::vector<float>& moments)
{
    assert(phases.size() == signal.size());

    moments.resize(n_moments + 1);
    compute_moments(
        phases.data(),
        phases.size(),
        signal.data(),
        n_moments,
        moments.data());
}


void solve_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& solution)
{
    solution.resize(first_column.size());
    solve_levinson(first_column.data(), solution.data(), first_column.size());
}


void dot_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& dot_product)
{
    dot_product.resize(first_column.size());
    dot_levinson(first_column.data(), dot_product.data(), first_column.size());
}


void levinson_from_dot(
    const std::vector<float>& dot_product,
    std::vector<float>& first_column)
{
    first_column.resize(dot_product.size());
    levinson_from_dot(dot_product.data(), first_column.data(), dot_product.size());
}


void compute_density(
    const std::vector<float>& phases, 
    const std::vector<float>& moments,
    std::vector<float>& density)
{
    density.resize(phases.size());
    compute_density(
        phases.data(), 
        phases.size(), 
        moments.data(), 
        moments.size() - 1, 
        density.data());
}


void compress_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments)
{
    dot_levinson(moments, compressed_moments);
    compressed_moments[0] = moments[0];
}


void decompress_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments)
{
    levinson_from_dot(compressed_moments, moments);
}

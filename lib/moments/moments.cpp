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

/******************************************************************************
 * C interface
 *****************************************************************************/

extern "C" {

void wavelengths_to_phases(
    const float wavelengths[],
    size_t n_wavelengths,
    float phases[])
{
    float min_wl = wavelengths[0];
    float max_wl = wavelengths[n_wavelengths - 1];

    for (size_t i = 0; i < n_wavelengths; i++) {
        phases[i] = M_PI * (wavelengths[i] - min_wl) / (max_wl - min_wl) - M_PI;
    }
}


void compute_moments(
    const float og_phases[],
    size_t n_phases,
    const float og_signal[],
    size_t n_moments, 
    float moments[])
{
    // TODO: double check, this looks a bit hacky
    std::vector<float> phases(n_phases + 2);
    std::vector<float> signal(n_phases + 2);

    memcpy(&phases[1], og_phases, n_phases * sizeof(float));
    memcpy(&signal[1], og_signal, n_phases * sizeof(float));

    phases[0] = -M_PI;
    signal[0] = og_signal[0];
    
    phases[n_phases] = 0;
    signal[n_phases] = og_signal[n_phases - 1];

    std::vector<std::complex<float>> t_moments(n_moments + 1);

    for (size_t i = 0; i < phases.size() - 1; i++) {
        // Required by the previous TODO
        if (phases[i] >= phases[i + 1]) {
            continue;
        }
 
        const float gradient = (signal[i + 1] - signal[i]) / (phases[i + 1] - phases[i]);
        const float y_intercept = signal[i] - gradient * phases[i];

        const std::complex<float> J = std::complex<float>(0., 1.);

        for (size_t k = 1; k < n_moments + 1; k++) {
            const std::complex<float> common_summands(
                gradient / (float)(k*k), 
                y_intercept / (float)k);

            t_moments[k] += 
                  (common_summands + gradient * J * (float)k * phases[i + 1] / (float)(k*k)) * std::exp(-J * (float)k * phases[i + 1])
                - (common_summands + gradient * J * (float)k * phases[i    ] / (float)(k*k)) * std::exp(-J * (float)k * phases[i    ]);
        }

        t_moments[0] += 
            (.5 * gradient * phases[i + 1] * phases[i + 1] + y_intercept * phases[i + 1])
          - (.5 * gradient * phases[i    ] * phases[i    ] + y_intercept * phases[i    ]);
     }

    // Mirrored signal - we keep the real part and multiply by 2.
    for (size_t k = 0; k < n_moments + 1; k++) {
        moments[k] = 1.f / M_PI * t_moments[k].real();
    }
}


void solve_levinson(
    const float first_column[],
    size_t size,
    float solution[])
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
    size_t size,
    float dot_product[])
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
    size_t size,
    float first_column[])
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
    dot_levinson(moments, n_moments + 1, compressed_moments);
    compressed_moments[0] = moments[0];
}


void decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[])
{
    levinson_from_dot(compressed_moments, n_compressed_moments + 1, moments);
}

} // extern "C"


/******************************************************************************
 * C++ interface
 *****************************************************************************/

void wavelengths_to_phases(
    const std::vector<float>& wavelengths, 
    std::vector<float>& phases) 
{
    phases.resize(wavelengths.size());
    wavelengths_to_phases(wavelengths.data(), wavelengths.size(), phases.data());
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
    solve_levinson(first_column.data(), first_column.size(), solution.data());
}


void dot_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& dot_product)
{
    dot_product.resize(first_column.size());
    dot_levinson(first_column.data(), first_column.size(), dot_product.data());
}


void levinson_from_dot(
    const std::vector<float>& dot_product,
    std::vector<float>& first_column)
{
    first_column.resize(dot_product.size());
    levinson_from_dot(dot_product.data(), dot_product.size(), first_column.data());
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

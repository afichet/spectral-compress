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

/******************************************************************************
 * C interface
 *****************************************************************************/

extern "C" {

void wavelengths_to_phases(
    const float wavelengths[],
    size_t n_wavelengths,
    float phases[])
{
    const float min_wl = wavelengths[0];
    const float max_wl = wavelengths[n_wavelengths - 1];

    for (size_t i = 0; i < n_wavelengths; i++) {
        phases[i] = M_PIf * (wavelengths[i] - min_wl) / (max_wl - min_wl) - M_PIf;
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

    // Cause a strange warning...
    // memcpy(&phases[1], og_phases, n_phases * sizeof(float));
    // memcpy(&signal[1], og_signal, n_phases * sizeof(float));

    for (size_t i = 0; i < n_phases; i++) {
        phases[i + 1] = og_phases[i];
        signal[i + 1] = og_signal[i];
    }

    phases[0]        = -M_PIf;
    phases[n_phases] = 0;

    signal[0]        = og_signal[0];
    signal[n_phases] = og_signal[n_phases - 1];

    std::vector<std::complex<float>> t_moments(n_moments + 1);
    const std::complex<float> J = std::complex<float>(0., 1.);

    for (size_t i = 0; i < phases.size() - 1; i++) {
        // Required by the previous TODO
        if (phases[i] >= phases[i + 1]) {
            continue;
        }
 
        const float gradient    = (signal[i + 1] - signal[i]) / (phases[i + 1] - phases[i]);
        const float y_intercept = signal[i] - gradient * phases[i];

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

    // Mirrored signal
    for (size_t k = 0; k < n_moments + 1; k++) {
        moments[k] = t_moments[k].real() / M_PIf;
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
        v[i] = moments[i] / (2.f * M_PIf);
    }

    std::vector<float> solution;
    solve_levinson(v, solution);

    std::vector<std::complex<float>> denum(n_phases);

    const float r = 1.f / (2.f * M_PIf);

    for (size_t p = 0; p < n_phases; p++) {
        std::complex<float> denum = 0;

        for (size_t k = 0; k < n_moments + 1; k++) {
            denum += solution[k] * std::polar(r, k * phases[p]);
        }

        denum = std::abs(denum);

        density[p] = r * solution[0] / (denum * denum).real();
    }
}


void compute_density_bounded_lagrange(
    const float phases[],
    size_t n_phases,
    const float moments[],
    size_t n_moments,
    float density[])
{
    std::vector<std::complex<float>> exponential_moments(n_moments + 1);
    std::vector<std::complex<float>> toeplitz_column(n_moments + 1);

    moments_to_exponential_moments(moments, n_moments, exponential_moments);

    for (size_t i = 1; i < n_moments + 1; i++) {
        toeplitz_column[i] = exponential_moments[i] / (2.f * M_PIf);
    }

    toeplitz_column[0] = exponential_moments[0].real() / M_PIf;

    std::vector<std::complex<float>> evaluation_polynomial;
    solve_levinson(toeplitz_column, evaluation_polynomial);

    std::vector<std::complex<float>> lagrange_multipliers;
    compute_lagrange_multipliers(exponential_moments, evaluation_polynomial, lagrange_multipliers);

    // Evaluate the Fourier series
    std::vector<float> fourier_series(n_phases);

    for (size_t i = 1; i < n_moments + 1; i++) {
        for (size_t p = 0; p < n_phases; p++) {
            fourier_series[p] += (
                lagrange_multipliers[i] * std::exp(std::complex<float>(0.f, (float)i * phases[p]))
                ).real();
        }
    }

    for (size_t p = 0; p < n_phases; p++) {
        fourier_series[p] = 2.f * fourier_series[p] + lagrange_multipliers[0].real();
        density[p] = std::atan(fourier_series[p]) / M_PIf + .5f;
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


void compress_bounded_moments(
    const float moments[],
    size_t n_moments,
    float compressed_moments[])
{
    std::vector<std::complex<float>> exponential_moments(n_moments + 1);

    moments_to_exponential_moments(
        moments,
        n_moments + 1,
        exponential_moments
    );

    std::vector<std::complex<float>> toeplitz_first_column(n_moments + 1);

    for (size_t i = 0; i < n_moments + 1; i++) {
        toeplitz_first_column[i] = exponential_moments[i] / (2.f * M_PIf);
    }

    toeplitz_first_column[0] = 2.f * toeplitz_first_column[0].real();

    std::vector<std::complex<float>> dots;
    dot_levinson(toeplitz_first_column, dots);

    const std::complex<float> m = std::abs(exponential_moments[0]) / (std::complex<float>(0.f, 1.f) * exponential_moments[0]);

    for (size_t i = 1; i < n_moments + 1; i++) {
        compressed_moments[i] = (dots[i] * m).real();
    }

    compressed_moments[0] = moments[0];
}


void decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[])
{
    levinson_from_dot(compressed_moments, n_compressed_moments + 1, moments);
}


void decompress_bounded_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[])
{
    const std::complex<float> J(0.f, 1.f);
    const std::complex<float> exp_0 = std::exp(J * M_PIf * (compressed_moments[0] - .5f)) / (4.f * M_PIf);

    std::vector<std::complex<float>> dots(n_compressed_moments + 1);

    for (size_t i = 1; i < n_compressed_moments + 1; i++) {
        dots[i] = compressed_moments[i] * J * exp_0 / std::abs(exp_0);
    }

    dots[0] = exp_0.real() / M_PIf;

    std::vector<std::complex<float>> exponential_moments;
    levinson_from_dot(dots, exponential_moments);

    for (size_t i = 1; i < n_compressed_moments + 1; i++) {
        exponential_moments[i] *= 2.f * M_PIf;
    }

    exponential_moments[0] = exp_0;

    exponential_moments_to_moments(exponential_moments, moments);
}


/* ----------------------------------------------------------------------------
   Utility functions
   ------------------------------------------------------------------------- */

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
        const float radius = 1.f / solution[0];
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


void compute_density_bounded_lagrange(
    const std::vector<float>& phases,
    const std::vector<float>& moments,
    std::vector<float>& density)
{
    density.resize(phases.size());
    compute_density_bounded_lagrange(
        phases.data(),
        phases.size(),
        moments.data(),
        moments.size() - 1,
        density.data()
    );
}


void compress_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments)
{
    dot_levinson(moments, compressed_moments);
    compressed_moments[0] = moments[0];
}


void compress_bounded_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments)
{
    compressed_moments.resize(moments.size());
    compress_bounded_moments(
        moments.data(),
        moments.size() - 1,
        compressed_moments.data()
    );
}


void decompress_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments)
{
    levinson_from_dot(compressed_moments, moments);
}


void decompress_bounded_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments)
{
    moments.resize(compressed_moments.size());
    decompress_bounded_moments(
        compressed_moments.data(),
        compressed_moments.size() - 1,
        moments.data()
    );
}


/* ----------------------------------------------------------------------------
   Utility functions
   ------------------------------------------------------------------------- */

void solve_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& solution)
{
    solution.resize(first_column.size());
    solve_levinson(first_column.data(), first_column.size(), solution.data());
}


void solve_levinson(
    const std::vector<std::complex<float>>& first_column,
    std::vector<std::complex<float>>& solution)
{
    solution.resize(first_column.size());
    solution[0] = 1.f / first_column[0];

    std::vector<std::complex<float>> temp_s(first_column.size());

    for (size_t i = 1; i < first_column.size(); i++) {
        solution[i] = 0;

        std::complex<float> dot_product = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product * std::conj(solution[i - k])) / (1.f - std::abs(dot_product) * std::abs(dot_product));
        }

        solution = temp_s;
    }
}


void dot_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& dot_product)
{
    dot_product.resize(first_column.size());
    dot_levinson(first_column.data(), first_column.size(), dot_product.data());
}


void dot_levinson(
    const std::vector<std::complex<float>>& first_column,
    std::vector<std::complex<float>>& dot_product)
{
    dot_product.resize(first_column.size());

    std::vector<std::complex<float>> solution(first_column.size());
    std::vector<std::complex<float>> temp_s(first_column.size());

    dot_product[0] = 0;
    solution[0]    = 1.f / first_column[0];

    for (size_t i = 1; i < first_column.size(); i++) {
        dot_product[i] = 0;
        solution[i]    = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product[i] += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * std::conj(solution[i - k])) / (1.f - std::abs(dot_product[i]) * std::abs(dot_product[i]));
        }

        solution = temp_s;
    }
}


void levinson_from_dot(
    const std::vector<float>& dot_product,
    std::vector<float>& first_column)
{
    first_column.resize(dot_product.size());
    levinson_from_dot(dot_product.data(), dot_product.size(), first_column.data());
}


void levinson_from_dot(
    const std::vector<std::complex<float>>& dot_product,
    std::vector<std::complex<float>>& first_column)
{
    first_column.resize(dot_product.size());

    std::vector<std::complex<float>> solution(dot_product.size());
    std::vector<std::complex<float>> temp_s(dot_product.size());

    first_column[0] = dot_product[0];
    solution[0]     = 1.f / dot_product[0];

    for (size_t i = 1; i < dot_product.size(); i++) {
        const float radius = 1.f / solution[0].real();
        std::complex<float> center = 0;

        for (size_t k = 1; k < i; k++) {
            center += solution[k] * first_column[i - k];
        }

        first_column[i] = radius * (dot_product[i] - center);

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * std::conj(solution[i - k])) / (1.f - std::abs(dot_product[i]) * std::abs(dot_product[i]));
        }

        solution = temp_s;
    }
}


void moments_to_exponential_moments(
    const float moments[],
    size_t n_moments,
    std::vector<std::complex<float>>& exponential_moments)
{
    exponential_moments.resize(n_moments + 1);

    exponential_moments[0] = std::exp(std::complex<float>(0.f, M_PIf * (moments[0] - .5f))) / (4.f * M_PIf);

    for (size_t i = 1; i < n_moments + 1; i++) {
        for (size_t k = 0; k < i; k++) {
            exponential_moments[i] += (float)(i - k) * exponential_moments[k] * moments[i - k];
        }

        exponential_moments[i] *= std::complex<float>(0.f, 2.f * M_PIf) / (float)i;
    }
}


void exponential_moments_to_moments(
    const std::vector<std::complex<float>>& exponential_moments,
    float moments[])
{
    const std::complex<float> J(0.f, 1.f);

    memset(moments, 0, exponential_moments.size() * sizeof(float));
    moments[0] = std::arg(exponential_moments[0]) / M_PIf + .5f;
    
    const std::complex<float> exp_0 = std::exp(J * M_PIf * (moments[0] - .5f)) / (4.f * M_PIf);

    for (size_t i = 1; i < exponential_moments.size(); i++) {
        std::complex<float> sum(0);

        for (size_t k = 1; k < i; k++) {
            sum += moments[i - k] * exponential_moments[k] * (float)(i - k);
        }
        
        moments[i] = (exponential_moments[i] / (J * 2.f * M_PIf * exp_0) - 1.f / ((float)i * exp_0) * sum).real();
    }
}


void compute_lagrange_multipliers(
    const std::vector<std::complex<float>>& exponential_moments,
    const std::vector<std::complex<float>>& evaluation_polynomial,
    std::vector<std::complex<float>>& lagrange_multipliers)
{
    std::vector<std::complex<float>> autocorrelation(evaluation_polynomial.size());

    for (size_t i = 0; i < exponential_moments.size(); i++) {
        autocorrelation[i] = 0;
        
        for (size_t k = 0; k < exponential_moments.size() - i; k++) {
            autocorrelation[i] += std::conj(evaluation_polynomial[i + k]) * evaluation_polynomial[k];
        }
    }

    lagrange_multipliers.resize(exponential_moments.size());

    for (size_t i = 0; i < exponential_moments.size(); i++) {
        lagrange_multipliers[i] = 0;

        for (size_t k = 0; k < exponential_moments.size() - i; k++) {
            lagrange_multipliers[i] += exponential_moments[k] * autocorrelation[k + i];
        }

        lagrange_multipliers[i] /= evaluation_polynomial[0] * std::complex<float>(0, M_PIf);
    }
}

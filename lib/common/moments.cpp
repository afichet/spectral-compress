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

#include "moments.h"

#include <cmath>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#define MIN_FLT 1e-35

/******************************************************************************
 * C interface
 *****************************************************************************/

extern "C"
void wavelengths_to_phases(
    const double wavelengths[],
    size_t n_wavelengths,
    double phases[])
{
    const double min_wl = wavelengths[0];
    const double max_wl = wavelengths[n_wavelengths - 1];

    for (size_t i = 0; i < n_wavelengths; i++) {
        phases[i] = M_PI * (wavelengths[i] - min_wl) / (max_wl - min_wl) - M_PI;
    }
}

void wavelengths_to_phases(
    const std::vector<double>& wavelengths,
    std::vector<double>& phases)
{
    phases.resize(wavelengths.size());
    wavelengths_to_phases(wavelengths.data(), wavelengths.size(), phases.data());
}


extern "C"
void compute_basis_signal_to_moments(
    const double phases[],
    size_t n_phases,
    double basis[])
{
    // Eigen::Map<Eigen::MatrixXf> transform_mat(basis, n_phases, n_phases);
    Eigen::MatrixXd signal(n_phases, n_phases);

    signal.setIdentity();

    #pragma omp parallel for
    for (size_t i = 0; i < n_phases; i++) {
        compute_moments(
            phases, n_phases,
            signal.col(i).data(),
            n_phases,
            &basis[i * n_phases]
        );
    }
}

void compute_basis_signal_to_moments(
    const std::vector<double>& phases,
    std::vector<double>& basis)
{
    basis.resize(phases.size() * phases.size());

    compute_basis_signal_to_moments(
        phases.data(), phases.size(),
        basis.data()
    );
}


extern "C"
void compute_basis_moments_to_signal(
    const double phases[],
    size_t n_phases,
    double basis[])
{
    compute_basis_signal_to_moments(phases, n_phases, basis);
    Eigen::Map<Eigen::MatrixXd> transform(basis, n_phases, n_phases);
    transform = transform.inverse();
}

void compute_basis_moments_to_signal(
    const std::vector<double>& phases,
    std::vector<double>& basis)
{
    basis.resize(phases.size() * phases.size());

    compute_basis_moments_to_signal(
        phases.data(), phases.size(),
        basis.data()
    );
}


extern "C"
void compute_moments(
    const double og_phases[],
    size_t n_phases,
    const double og_signal[],
    size_t n_moments,
    double moments[])
{
    assert(og_phases[0] >= -M_PI);
    assert(og_phases[n_moments - 1] <= 0);

    // Duplicate first and last element if phases do not cover the whole
    // integration domain
    size_t vec_sz = n_phases;
    size_t start_idx = 0;

    if (og_phases[0] != -M_PI) {
        vec_sz++;
        start_idx = 1;
    }

    if (og_phases[n_moments - 1] < 0) {
        vec_sz++;
    }

    std::vector<double> phases(vec_sz);
    std::vector<double> signal(vec_sz);

    phases.front() = -M_PI;
    phases.back()  = 0;

    signal.front() = og_signal[0];
    signal.back()  = og_signal[n_phases - 1];

    for (size_t i = start_idx; i < n_phases; i++) {
        phases[i] = og_phases[i - start_idx];
        signal[i] = og_signal[i - start_idx];

        assert(phases[i] >= -M_PI);
        assert(signal[i] >= 0);
    }

    std::vector<std::complex<double>> t_moments(n_moments);
    const std::complex<double> J = std::complex<double>(0., 1.);

    for (size_t i = 0; i < phases.size() - 1; i++) {
        assert(phases[i] != phases[i + 1]);

        const double d_signal = signal[i + 1] - signal[i];
        const double d_phases = phases[i + 1] - phases[i];
        const double gradient = d_signal / d_phases;
        const double y_intercept = signal[i] - gradient * phases[i];

        for (size_t k = 1; k < n_moments; k++) {
            const double k_d = (double)k;
            // const std::complex<double> common_summands(
            //     gradient / (double)(k*k),
            //     y_intercept / (double)k);

            const std::complex<double> common_summands =
                std::complex<double>(gradient / k_d, y_intercept) / k_d;

            t_moments[k] +=
                  (common_summands + gradient * J * phases[i + 1] / k_d) * std::exp(-J * k_d * phases[i + 1])
                - (common_summands + gradient * J * phases[i    ] / k_d) * std::exp(-J * k_d * phases[i    ]);
        }

        // t_moments[0] +=
        //     (.5 * gradient * phases[i + 1] * phases[i + 1] + y_intercept * phases[i + 1])
        //   - (.5 * gradient * phases[i    ] * phases[i    ] + y_intercept * phases[i    ]);

        // t_moments[0] += .5 * d_signal * (phases[i + 1] + phases[i]) + y_intercept * d_phases;
        t_moments[0] +=
            d_signal * (phases[i + 1] + phases[i]) / 2.
            + signal[i] * d_phases
            - d_signal * phases[i];
     }

    // Mirrored signal
    for (size_t k = 0; k < n_moments; k++) {
        moments[k] = t_moments[k].real() / M_PI;
    }
}

void compute_moments(
    const std::vector<double>& phases,
    const std::vector<double>& signal,
    size_t n_moments,
    std::vector<double>& moments)
{
    assert(phases.size() == signal.size());

    moments.resize(n_moments);

    compute_moments(
        phases.data(),
        phases.size(),
        signal.data(),
        n_moments,
        moments.data());
}


extern "C"
void compute_density(
    const double phases[],
    size_t n_phases,
    const double moments[],
    size_t n_moments,
    double density[])
{
    if (moments[0] > 0) {
        std::vector<double> v(n_moments);

        for (size_t i = 0; i < v.size(); i++) {
            v[i] = moments[i] / (2. * M_PI);
        }

        std::vector<double> solution;
        solve_levinson(v, solution);

        std::vector<std::complex<double>> denum(n_phases);

        const double r = 1. / (2. * M_PI);

        for (size_t p = 0; p < n_phases; p++) {
            std::complex<double> denum = 0;

            for (size_t k = 0; k < n_moments; k++) {
                denum += solution[k] * std::polar(r, k * phases[p]);
            }

            denum = std::abs(denum);

            density[p] = r * solution[0] / (denum * denum).real();
        }
    } else {
        for (size_t p = 0; p < n_phases; p++) {
            density[p] = 0;
        }
    }
}

void compute_density(
    const std::vector<double>& phases,
    const std::vector<double>& moments,
    std::vector<double>& density)
{
    density.resize(phases.size());

    compute_density(
        phases.data(),
        phases.size(),
        moments.data(),
        moments.size(),
        density.data());
}


extern "C"
void bounded_compute_density_lagrange(
    const double phases[],
    size_t n_phases,
    const double moments[],
    size_t n_moments,
    double density[])
{
    if (moments[0] > 0) {
        std::vector<std::complex<double>> exponential_moments(n_moments);
        std::vector<std::complex<double>> toeplitz_column(n_moments);

        moments_to_exponential_moments(moments, n_moments, exponential_moments);

        for (size_t i = 1; i < n_moments; i++) {
            toeplitz_column[i] = exponential_moments[i] / (2. * M_PI);
        }

        toeplitz_column[0] = exponential_moments[0].real() / M_PI;

        std::vector<std::complex<double>> evaluation_polynomial;
        solve_levinson(toeplitz_column, evaluation_polynomial);

        std::vector<std::complex<double>> lagrange_multipliers;
        compute_lagrange_multipliers(exponential_moments, evaluation_polynomial, lagrange_multipliers);

        // Evaluate the Fourier series
        std::vector<double> fourier_series(n_phases);

        for (size_t i = 1; i < n_moments; i++) {
            for (size_t p = 0; p < n_phases; p++) {
                fourier_series[p] += (
                    lagrange_multipliers[i] * std::exp(std::complex<double>(0., (double)i * phases[p]))
                    ).real();
            }
        }

        for (size_t p = 0; p < n_phases; p++) {
            fourier_series[p] = 2. * fourier_series[p] + lagrange_multipliers[0].real();
            density[p] = std::atan(fourier_series[p]) / M_PI + .5f;
        }
    } else {
        for (size_t p = 0; p < n_phases; p++) {
            density[p] = 0;
        }
    }
}

void bounded_compute_density_lagrange(
    const std::vector<double>& phases,
    const std::vector<double>& moments,
    std::vector<double>& density)
{
    density.resize(phases.size());

    bounded_compute_density_lagrange(
        phases.data(),
        phases.size(),
        moments.data(),
        moments.size(),
        density.data()
    );
}


/*****************************************************************************/
/* Compression                                                               */
/*****************************************************************************/

extern "C"
void unbounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[])
{
    if (moments[0] > 0) {
        dot_levinson(moments, n_moments, compressed_moments);
        compressed_moments[0] = moments[0];
    } else {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments[m] = 0;
        }
    }
}

void unbounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments)
{
    dot_levinson(moments, compressed_moments);
    compressed_moments[0] = moments[0];
}


extern "C"
void bounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[])
{
    if (moments[0] > 0) {
        std::vector<std::complex<double>> exponential_moments(n_moments);

        moments_to_exponential_moments(
            moments,
            n_moments,
            exponential_moments
        );

        std::vector<std::complex<double>> toeplitz_first_column(n_moments);

        for (size_t i = 0; i < n_moments; i++) {
            toeplitz_first_column[i] = exponential_moments[i] / (2. * M_PI);
        }

        toeplitz_first_column[0] = 2. * toeplitz_first_column[0].real();

        std::vector<std::complex<double>> dots;
        dot_levinson(toeplitz_first_column, dots);

        const std::complex<double> m = std::abs(exponential_moments[0]) / (std::complex<double>(0., 1.) * exponential_moments[0]);
        assert(!std::isinf(m.real()));
        assert(!std::isinf(m.imag()));
        assert(!std::isnan(m.real()));
        assert(!std::isnan(m.imag()));

        for (size_t i = 1; i < n_moments; i++) {
            compressed_moments[i] = (dots[i] * m).real();
            assert(!std::isinf(compressed_moments[i]));
            assert(!std::isnan(compressed_moments[i]));
        }

        compressed_moments[0] = moments[0];
        assert(!std::isinf(compressed_moments[0]));
        assert(!std::isnan(compressed_moments[0]));
    } else {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments[m] = 0;
        }
    }
}

void bounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments)
{
    compressed_moments.resize(moments.size());

    bounded_compress_moments(
        moments.data(),
        moments.size(),
        compressed_moments.data()
    );
}


extern "C"
void unbounded_to_bounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[])
{
    if (moments[0] > 0) {
        std::vector<double> rescaled_moments(n_moments);
        std::vector<std::complex<double>> exponential_moments(n_moments);

        for (size_t i = 0; i < n_moments; i++) {
            rescaled_moments[i] = moments[i] / ((double)n_moments * moments[0]);
        }

        moments_to_exponential_moments(
            rescaled_moments.data(),
            n_moments,
            exponential_moments
        );

        std::vector<std::complex<double>> toeplitz_first_column(n_moments);

        for (size_t i = 0; i < n_moments; i++) {
            toeplitz_first_column[i] = exponential_moments[i] / (2. * M_PI);
        }

        toeplitz_first_column[0] = 2. * toeplitz_first_column[0].real();

        std::vector<std::complex<double>> dots;
        dot_levinson(toeplitz_first_column, dots);

        const std::complex<double> m = std::abs(exponential_moments[0]) / (std::complex<double>(0., 1.) * exponential_moments[0]);

        for (size_t i = 1; i < n_moments; i++) {
            compressed_moments[i] = (dots[i] * m).real();
        }

        compressed_moments[0] = moments[0];
    } else {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments[m] = 0;
        }
    }
}

void unbounded_to_bounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments)
{
    compressed_moments.resize(moments.size());

    unbounded_to_bounded_compress_moments(
        moments.data(),
        moments.size(),
        compressed_moments.data()
    );
}


/*****************************************************************************/
/* Decompression                                                             */
/*****************************************************************************/

extern "C"
void unbounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[])
{
    if (compressed_moments[0] > 0) {
        levinson_from_dot(compressed_moments, n_compressed_moments, moments);
    } else {
        for (size_t m = 0; m < n_compressed_moments; m++) {
            moments[m] = 0;
        }
    }
}

void unbounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments)
{
    levinson_from_dot(compressed_moments, moments);
}



extern "C"
void bounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[])
{
    if (compressed_moments[0] > 0) {
        const std::complex<double> J(0., 1.);
        const std::complex<double> exp_0 = std::exp(J * M_PI * (compressed_moments[0] - .5)) / (4. * M_PI);
        assert(!std::isinf(exp_0.real()));
        assert(!std::isinf(exp_0.imag()));
        assert(!std::isnan(exp_0.real()));
        assert(!std::isnan(exp_0.imag()));

        std::vector<std::complex<double>> dots(n_compressed_moments);

        for (size_t i = 1; i < n_compressed_moments; i++) {
            dots[i] = compressed_moments[i] * J * exp_0 / std::abs(exp_0);
            assert(!std::isinf(dots[i].real()));
            assert(!std::isinf(dots[i].imag()));
            assert(!std::isnan(dots[i].real()));
            assert(!std::isnan(dots[i].imag()));
        }

        dots[0] = exp_0.real() / M_PI;

        std::vector<std::complex<double>> exponential_moments;
        levinson_from_dot(dots, exponential_moments);

        for (size_t i = 1; i < n_compressed_moments; i++) {
            exponential_moments[i] *= 2. * M_PI;
        }

        exponential_moments[0] = exp_0;

        exponential_moments_to_moments(exponential_moments, moments);
    } else {
        for (size_t m = 0; m < n_compressed_moments; m++) {
            moments[m] = 0;
        }
    }
}

void bounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments)
{
    moments.resize(compressed_moments.size());
    bounded_decompress_moments(
        compressed_moments.data(),
        compressed_moments.size(),
        moments.data()
    );
}


extern "C"
void unbounded_to_bounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[])
{
    if (compressed_moments[0] > 0) {
        const std::complex<double> J(0., 1.);
        const double m0 = 1. / (double)n_compressed_moments;

        const std::complex<double> exp_0 = std::exp(J * M_PI * (m0 - .5)) / (4. * M_PI);

        std::vector<std::complex<double>> dots(n_compressed_moments);

        for (size_t i = 1; i < n_compressed_moments; i++) {
            dots[i] = compressed_moments[i] * J * exp_0 / std::abs(exp_0);
        }

        dots[0] = exp_0.real() / M_PI;

        std::vector<std::complex<double>> exponential_moments;
        levinson_from_dot(dots, exponential_moments);

        for (size_t i = 1; i < n_compressed_moments; i++) {
            exponential_moments[i] *= 2. * M_PI;
        }

        exponential_moments[0] = exp_0;

        exponential_moments_to_moments(exponential_moments, moments);

        for (size_t i = 0; i < n_compressed_moments; i++) {
            moments[i] *= (double)n_compressed_moments * compressed_moments[0];
        }
    } else {
        for (size_t m = 0; m < n_compressed_moments; m++) {
            moments[m] = 0;
        }
    }
}

void unbounded_to_bounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments)
{
    moments.resize(compressed_moments.size());

    unbounded_to_bounded_decompress_moments(
        compressed_moments.data(),
        compressed_moments.size(),
        moments.data()
    );
}


/*****************************************************************************/
/* Utility functions                                                         */
/*****************************************************************************/

extern "C"
void solve_levinson(
    const double first_column[],
    size_t size,
    double solution[])
{
    solution[0] = 1. / first_column[0];

    std::vector<double> temp_s(size);

    for (size_t i = 1; i < size; i++) {
        solution[i] = 0;

        double dot_product = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product * solution[i - k]) / (1. - dot_product * dot_product);
        }

        memcpy(solution, temp_s.data(), temp_s.size() * sizeof(double));
    }
}

void solve_levinson(
    const std::vector<double>& first_column,
    std::vector<double>& solution)
{
    solution.resize(first_column.size());
    solve_levinson(first_column.data(), first_column.size(), solution.data());
}


extern "C"
void dot_levinson(
    const double first_column[],
    size_t size,
    double dot_product[])
{
    std::vector<double> solution(size);
    std::vector<double> temp_s(size);

    dot_product[0] = 0;
    solution[0]    = 1. / first_column[0];

    for (size_t i = 1; i < size; i++) {
        dot_product[i] = 0;
        solution[i]    = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product[i] += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * solution[i - k]) / (1. - dot_product[i] * dot_product[i]);
        }

        solution = temp_s;
    }
}

void dot_levinson(
    const std::vector<double>& first_column,
    std::vector<double>& dot_product)
{
    dot_product.resize(first_column.size());
    dot_levinson(first_column.data(), first_column.size(), dot_product.data());
}



extern "C"
void levinson_from_dot(
    const double dot_product[],
    size_t size,
    double first_column[])
{
    std::vector<double> solution(size);
    std::vector<double> temp_s(size);

    first_column[0] = dot_product[0];
    solution[0]     = 1. / dot_product[0];

    for (size_t i = 1; i < size; i++) {
        const double radius = 1. / solution[0];
        double center = 0;

        for (size_t k = 1; k < i; k++) {
            center += solution[k] * first_column[i - k];
        }

        first_column[i] = radius * (dot_product[i] - center);

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * solution[i - k]) / (1. - dot_product[i] * dot_product[i]);
        }

        solution = temp_s;
    }
}

void levinson_from_dot(
    const std::vector<double>& dot_product,
    std::vector<double>& first_column)
{
    first_column.resize(dot_product.size());
    levinson_from_dot(dot_product.data(), dot_product.size(), first_column.data());
}


void solve_levinson(
    const std::vector<std::complex<double>>& first_column,
    std::vector<std::complex<double>>& solution)
{
    solution.resize(first_column.size());
    solution[0] = 1. / first_column[0];

    std::vector<std::complex<double>> temp_s(first_column.size());

    for (size_t i = 1; i < first_column.size(); i++) {
        solution[i] = 0;

        std::complex<double> dot_product = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product * std::conj(solution[i - k])) / (1. - std::abs(dot_product) * std::abs(dot_product));
        }

        solution = temp_s;
    }
}


void dot_levinson(
    const std::vector<std::complex<double>>& first_column,
    std::vector<std::complex<double>>& dot_product)
{
    dot_product.resize(first_column.size());

    std::vector<std::complex<double>> solution(first_column.size());
    std::vector<std::complex<double>> temp_s(first_column.size());

    dot_product[0] = 0;
    solution[0]    = 1. / first_column[0];

    for (size_t i = 1; i < first_column.size(); i++) {
        dot_product[i] = 0;
        solution[i]    = 0;

        for (size_t k = 0; k < i; k++) {
            dot_product[i] += solution[k] * first_column[i - k];
        }

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * std::conj(solution[i - k])) / (1. - std::abs(dot_product[i]) * std::abs(dot_product[i]));
        }

        solution = temp_s;
    }
}


void levinson_from_dot(
    const std::vector<std::complex<double>>& dot_product,
    std::vector<std::complex<double>>& first_column)
{
    first_column.resize(dot_product.size());

    std::vector<std::complex<double>> solution(dot_product.size());
    std::vector<std::complex<double>> temp_s(dot_product.size());

    first_column[0] = dot_product[0];
    solution[0]     = 1. / dot_product[0];

    for (size_t i = 1; i < dot_product.size(); i++) {
        const double radius = 1. / solution[0].real();
        std::complex<double> center = 0;

        for (size_t k = 1; k < i; k++) {
            center += solution[k] * first_column[i - k];
        }

        first_column[i] = radius * (dot_product[i] - center);

        for (size_t k = 0; k < i + 1; k++) {
            temp_s[k] = (solution[k] - dot_product[i] * std::conj(solution[i - k])) / (1. - std::abs(dot_product[i]) * std::abs(dot_product[i]));
        }

        solution = temp_s;
    }
}


void moments_to_exponential_moments(
    const double moments[],
    size_t n_moments,
    std::vector<std::complex<double>>& exponential_moments)
{
    exponential_moments.resize(n_moments);

    exponential_moments[0] = std::exp(std::complex<double>(0., M_PI * (moments[0] - .5))) / (4. * M_PI);
    assert(!std::isnan(exponential_moments[0].real()));
    assert(!std::isnan(exponential_moments[0].imag()));
    assert(!std::isinf(exponential_moments[0].real()));
    assert(!std::isinf(exponential_moments[0].imag()));

    for (size_t i = 1; i < n_moments; i++) {
        for (size_t k = 0; k < i; k++) {
            exponential_moments[i] += (double)(i - k) * exponential_moments[k] * moments[i - k];
        }

        exponential_moments[i] *= std::complex<double>(0., 2. * M_PI) / (double)i;
        assert(!std::isnan(exponential_moments[i].real()));
        assert(!std::isnan(exponential_moments[i].imag()));
        assert(!std::isinf(exponential_moments[i].real()));
        assert(!std::isinf(exponential_moments[i].imag()));
    }
}


void exponential_moments_to_moments(
    const std::vector<std::complex<double>>& exponential_moments,
    double moments[])
{
    const std::complex<double> J(0., 1.);

    memset(moments, 0, exponential_moments.size() * sizeof(double));
    moments[0] = std::arg(exponential_moments[0]) / M_PI + .5;

    const std::complex<double> exp_0 = std::exp(J * M_PI * (moments[0] - .5)) / (4. * M_PI);

    for (size_t i = 1; i < exponential_moments.size(); i++) {
        std::complex<double> sum(0);

        for (size_t k = 1; k < i; k++) {
            sum += moments[i - k] * exponential_moments[k] * (double)(i - k);
        }

        moments[i] = (exponential_moments[i] / (J * 2. * M_PI * exp_0) - 1. / ((double)i * exp_0) * sum).real();
    }
}


void compute_lagrange_multipliers(
    const std::vector<std::complex<double>>& exponential_moments,
    const std::vector<std::complex<double>>& evaluation_polynomial,
    std::vector<std::complex<double>>& lagrange_multipliers)
{
    std::vector<std::complex<double>> autocorrelation(evaluation_polynomial.size());

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

        lagrange_multipliers[i] /= evaluation_polynomial[0] * std::complex<double>(0, M_PI);
    }
}

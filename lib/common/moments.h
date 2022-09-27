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

/******************************************************************************
 * C interface
 *****************************************************************************/

extern "C" {

#include <stddef.h>

/**
 * @brief Computes corresponding phases from given wavelengths
 * 
 * @param wavelengths   Wavelenghts to be converted to phases.
 * @param n_wavelengths Number of wavelengths in the array.
 * @param phases        Computed phases. Must be allocated with `n_wavelength` 
 *                      elements.
 */
void wavelengths_to_phases(
    const float wavelengths[],
    size_t n_wavelengths,
    float phases[]);


void compute_basis_signal_to_moments(
    const float phases[],
    size_t n_phases,
    float basis[]);


void compute_basis_moments_to_signal(
    const float phases[],
    size_t n_phases,
    float basis[]);


/**
 * @brief Computes moments from a given signal
 * 
 * @param phases    Phases at which the signal was sampled.
 * @param n_phases  Number of samples.
 * @param signal    Sampled value of the signal at the corresponding phases.
 * @param n_moments Number of moments to compute.
 * @param moments   Computed array of moments. Must be allocated with
 *                  `n_moments` elements.
 */
void compute_moments(
    const float phases[],
    size_t n_phases,
    const float signal[],
    size_t n_moments, 
    float moments[]);


/**
 * @brief Computes a density matching the given moments
 * 
 * @param phases    Phases where the density shall be computed.
 * @param n_phases  Size of phases array.
 * @param moments   Trigonometric moments.
 * @param n_moments Number of trigonometric moments.
 * @param density   Computed density. Must be allocated with `n_phases` 
 *                  elements.
 */
void compute_density(
    const float phases[],
    size_t n_phases,
    const float moments[],
    size_t n_moments,
    float density[]);


/**
 * @brief Computes a bounded density matching the given moments
 * 
 * @param phases    Phases where the density shall be computed.
 * @param n_phases  Size of phases array.
 * @param moments   Trigonometric moments.
 * @param n_moments Number of trigonometric moments.
 * @param density   Computed density. Must be allocated with `n_phases` 
 *                  elements.
 */
void bounded_compute_density_lagrange(
    const float phases[],
    size_t n_phases,
    const float moments[],
    size_t n_moments,
    float density[]);


/**
 * @brief Compress a set of moment
 * 
 * @param moments             Moments to compress.
 * @param n_moments           Number of moments.
 * @param compressed_moments  Computed compressed moments. Must be allocated
 *                            with `n_compressed_moments` elements.
 */
void unbounded_compress_moments(
    const float moments[],
    size_t n_moments,
    float compressed_moments[]);


/**
 * @brief Compress a set of bounded moment
 * 
 * @param moments             Moments to compress.
 * @param n_moments           Number of moments.
 * @param compressed_moments  Computed compressed moments. Must be allocated
 *                            with `n_compressed_moments` elements.
 */
void bounded_compress_moments(
    const float moments[],
    size_t n_moments,
    float compressed_moments[]);


/**
 * @brief Decompress a set of compressed moments.
 * 
 * @param compressed_moments   Compressed moments.
 * @param n_compressed_moments Number of compressed moments.
 * @param moments              Computed moments. Must be allocated with 
 *                             `n_compressed_moments` elements.
 */
void unbounded_decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[]);


/**
 * @brief Decompress a set of compressed bounded moments.
 * 
 * @param compressed_moments   Compressed moments.
 * @param n_compressed_moments Number of compressed moments.
 * @param moments              Computed moments. Must be allocated with 
 *                             `n_compressed_moments` elements.
 */
void bounded_decompress_moments(
    const float compressed_moments[],
    size_t n_compressed_moments,
    float moments[]);


/* ----------------------------------------------------------------------------
   Utility functions
   ------------------------------------------------------------------------- */

/**
 * @brief Use the Levinson algorithm to solve Toepliz matrix multiplied with
 * unit vector
 * 
 * @param first_column First column of the Toeplitz matrix.
 * @param size         Size of the column.
 * @param solution     Solution of the equation. Must be allocated with `size`
 *                     elements.
 */
void solve_levinson(
    const float first_column[],
    size_t size,
    float solution[]);


/**
 * @brief Get the dot products computed from the Levinson algorithm
 * 
 * @param first_column First column of the Toeplitz matrix.
 * @param size         Size of the column.
 * @param dot_product  Dot product computed by the Levinson algorithm. Must be
 *                     allocated with `size` elements.
 */
void dot_levinson(
    const float first_column[],
    size_t size,
    float dot_product[]);


/**
 * @brief Computes the solution given the dot products
 * 
 * @param dot_product  Dot products.
 * @param size         Size of the dot products.
 * @param first_column Solution. Must be allocated with `size` elements.
 */
void levinson_from_dot(
    const float dot_product[],
    size_t size,
    float first_column[]);

} // extern "C"


/******************************************************************************
 * C++ interface
 *****************************************************************************/

#ifdef __cplusplus

#include <vector>
#include <complex>


void wavelengths_to_phases(
    const std::vector<float>& wavelengths, 
    std::vector<float>& phases);


void compute_basis_signal_to_moments(
    const std::vector<float>& phases,
    std::vector<float>& basis);


void compute_basis_moments_to_signal(
    const std::vector<float>& phases,
    std::vector<float>& basis);


void compute_moments(
    const std::vector<float>& phases, 
    const std::vector<float>& signal, 
    size_t n_moments,
    std::vector<float>& moments);


void compute_density(
    const std::vector<float>& phases, 
    const std::vector<float>& moments,
    std::vector<float>& density);


void bounded_compute_density_lagrange(
    const std::vector<float>& phases,
    const std::vector<float>& moments,
    std::vector<float>& density);


void unbounded_compress_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments);


void bounded_compress_moments(
    const std::vector<float>& moments,
    std::vector<float>& compressed_moments);


void unbounded_decompress_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments);


void bounded_decompress_moments(
    const std::vector<float>& compressed_moments,
    std::vector<float>& moments);

// void bounded_compress_moments(
//     const std::vector<float>& moments,
//     std::vector<float>&c compressed_moments);


/* ----------------------------------------------------------------------------
   Utility functions
   ------------------------------------------------------------------------- */

void solve_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& solution);


void solve_levinson(
    const std::vector<std::complex<float>>& first_column,
    std::vector<std::complex<float>>& solution);


void dot_levinson(
    const std::vector<float>& first_column,
    std::vector<float>& dot_product);


void dot_levinson(
    const std::vector<std::complex<float>>& first_column,
    std::vector<std::complex<float>>& dot_product);


void levinson_from_dot(
    const std::vector<float>& dot_product,
    std::vector<float>& first_column);


void levinson_from_dot(
    const std::vector<std::complex<float>>& dot_product,
    std::vector<std::complex<float>>& first_column);


void moments_to_exponential_moments(
    const float moments[],
    size_t n_moments,
    std::vector<std::complex<float>>& exponential_moments);


void exponential_moments_to_moments(
    const std::vector<std::complex<float>>& exponential_moments,
    float moments[]);


void compute_lagrange_multipliers(
    const std::vector<std::complex<float>>& exponential_moments,
    const std::vector<std::complex<float>>& evaluation_polynomial,
    std::vector<std::complex<float>>& lagrange_multipliers);



// void evaluate_herglotz_transform(
//     const std::vector<std::complex<float>>& point_in_disk,
//     const std::vector<std::complex<float>>& exponential_moments,
//     std::vector<std::complex<float>>& evaluation_polynomial)
// {
//     std::vector<std::complex<float>> coefficients(exponential_moments.size());

//     for (size_t i = exponential_moments.size() - 1; i >= 0; i--) {
//         coefficients[i] = evaluation_polynomial[exponential_moments.size() - 1 - i]
//             + coefficients[i + 1] / 
//     }
// }

#endif // __cplusplus



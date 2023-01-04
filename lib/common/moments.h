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

#pragma once

#include <stddef.h>

#ifdef __cplusplus
#include <vector>
#include <complex>
#endif // __cplusplus


/**
 * @brief Computes corresponding phases from given wavelengths
 *
 * @param wavelengths   Wavelenghts to be converted to phases.
 * @param n_wavelengths Number of wavelengths in the array.
 * @param phases        Computed phases. Must be allocated with `n_wavelength`
 *                      elements.
 */
extern "C"
void wavelengths_to_phases(
    const double wavelengths[],
    size_t n_wavelengths,
    double phases[]);

#ifdef __cplusplus
void wavelengths_to_phases(
    const std::vector<double>& wavelengths,
    std::vector<double>& phases);
#endif // __cplusplus

extern "C"
void compute_basis_signal_to_moments(
    const double phases[],
    size_t n_phases,
    double basis[]);

#ifdef __cplusplus
void compute_basis_signal_to_moments(
    const std::vector<double>& phases,
    std::vector<double>& basis);
#endif // __cplusplus

extern "C"
void compute_basis_moments_to_signal(
    const double phases[],
    size_t n_phases,
    double basis[]);

#ifdef __cplusplus
void compute_basis_moments_to_signal(
    const std::vector<double>& phases,
    std::vector<double>& basis);
#endif // __cplusplus

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
extern "C"
void compute_moments(
    const double phases[],
    size_t n_phases,
    const double signal[],
    size_t n_moments,
    double moments[]);

#ifdef __cplusplus
void compute_moments(
    const std::vector<double>& phases,
    const std::vector<double>& signal,
    size_t n_moments,
    std::vector<double>& moments);
#endif // __cplusplus

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
extern "C"
void compute_density(
    const double phases[],
    size_t n_phases,
    const double moments[],
    size_t n_moments,
    double density[]);

#ifdef __cplusplus
void compute_density(
    const std::vector<double>& phases,
    const std::vector<double>& moments,
    std::vector<double>& density);
#endif // __cplusplus

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
extern "C"
void bounded_compute_density_lagrange(
    const double phases[],
    size_t n_phases,
    const double moments[],
    size_t n_moments,
    double density[]);

#ifdef __cplusplus
void bounded_compute_density_lagrange(
    const std::vector<double>& phases,
    const std::vector<double>& moments,
    std::vector<double>& density);
#endif // __cplusplus

/*****************************************************************************/
/* Compression                                                               */
/*****************************************************************************/

/**
 * @brief Compress a set of moment
 *
 * @param moments             Moments to compress.
 * @param n_moments           Number of moments.
 * @param compressed_moments  Computed compressed moments. Must be allocated
 *                            with `n_compressed_moments` elements.
 */
extern "C"
void unbounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[]);

#ifdef __cplusplus
void unbounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments);
#endif // __cplusplus

/**
 * @brief Compress a set of bounded moment
 *
 * @param moments             Moments to compress.
 * @param n_moments           Number of moments.
 * @param compressed_moments  Computed compressed moments. Must be allocated
 *                            with `n_compressed_moments` elements.
 */
extern "C"
void bounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[]);

#ifdef __cplusplus
void bounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments);
#endif // __cplusplus

/**
 * @brief Compress a set of moment
 *
 * This uses zeroth moment to rescale the other moment in order to
 * use the bounded compression method
 *
 * @param moments             Moments to compress.
 * @param n_moments           Number of moments.
 * @param compressed_moments  Computed compressed moments. Must be allocated
 *                            with `n_compressed_moments` elements.
 */
extern "C"
void unbounded_to_bounded_compress_moments(
    const double moments[],
    size_t n_moments,
    double compressed_moments[]);

#ifdef __cplusplus
void unbounded_to_bounded_compress_moments(
    const std::vector<double>& moments,
    std::vector<double>& compressed_moments);
#endif // __cplusplus

/*****************************************************************************/
/* Decompression                                                             */
/*****************************************************************************/

/**
 * @brief Decompress a set of compressed moments.
 *
 * @param compressed_moments   Compressed moments.
 * @param n_compressed_moments Number of compressed moments.
 * @param moments              Computed moments. Must be allocated with
 *                             `n_compressed_moments` elements.
 */
extern "C"
void unbounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[]);

#ifdef __cplusplus
void unbounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments);
#endif // __cplusplus


/**
 * @brief Decompress a set of compressed bounded moments.
 *
 * @param compressed_moments   Compressed moments.
 * @param n_compressed_moments Number of compressed moments.
 * @param moments              Computed moments. Must be allocated with
 *                             `n_compressed_moments` elements.
 */
extern "C"
void bounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[]);

#ifdef __cplusplus
void bounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments);
#endif // __cplusplus


/**
 * @brief Decompress a set of compressed moments.
 *
 * @param compressed_moments   Compressed moments.
 * @param n_compressed_moments Number of compressed moments.
 * @param moments              Computed moments. Must be allocated with
 *                             `n_compressed_moments` elements.
 */
extern "C"
void unbounded_to_bounded_decompress_moments(
    const double compressed_moments[],
    size_t n_compressed_moments,
    double moments[]);

#ifdef __cplusplus
void unbounded_to_bounded_decompress_moments(
    const std::vector<double>& compressed_moments,
    std::vector<double>& moments);
#endif // __cplusplus

/*****************************************************************************/


/*****************************************************************************/
/* Utility functions                                                         */
/*****************************************************************************/

/**
 * @brief Use the Levinson algorithm to solve Toepliz matrix multiplied with
 * unit vector
 *
 * @param first_column First column of the Toeplitz matrix.
 * @param size         Size of the column.
 * @param solution     Solution of the equation. Must be allocated with `size`
 *                     elements.
 */
extern "C"
void solve_levinson(
    const double first_column[],
    size_t size,
    double solution[]);

#ifdef __cplusplus
void solve_levinson(
    const std::vector<double>& first_column,
    std::vector<double>& solution);
#endif // __cplusplus

/**
 * @brief Get the dot products computed from the Levinson algorithm
 *
 * @param first_column First column of the Toeplitz matrix.
 * @param size         Size of the column.
 * @param dot_product  Dot product computed by the Levinson algorithm. Must be
 *                     allocated with `size` elements.
 */
extern "C"
void dot_levinson(
    const double first_column[],
    size_t size,
    double dot_product[]);

#ifdef __cplusplus
void solve_levinson(
    const std::vector<std::complex<double>>& first_column,
    std::vector<std::complex<double>>& solution);
#endif // __cplusplus

/**
 * @brief Computes the solution given the dot products
 *
 * @param dot_product  Dot products.
 * @param size         Size of the dot products.
 * @param first_column Solution. Must be allocated with `size` elements.
 */
extern "C"
void levinson_from_dot(
    const double dot_product[],
    size_t size,
    double first_column[]);

#ifdef __cplusplus
void levinson_from_dot(
    const std::vector<double>& dot_product,
    std::vector<double>& first_column);
#endif // __cplusplus




#ifdef __cplusplus
void dot_levinson(
    const std::vector<double>& first_column,
    std::vector<double>& dot_product);

void dot_levinson(
    const std::vector<std::complex<double>>& first_column,
    std::vector<std::complex<double>>& dot_product);

void levinson_from_dot(
    const std::vector<std::complex<double>>& dot_product,
    std::vector<std::complex<double>>& first_column);

void moments_to_exponential_moments(
    const double moments[],
    size_t n_moments,
    std::vector<std::complex<double>>& exponential_moments);

void exponential_moments_to_moments(
    const std::vector<std::complex<double>>& exponential_moments,
    double moments[]);

void compute_lagrange_multipliers(
    const std::vector<std::complex<double>>& exponential_moments,
    const std::vector<std::complex<double>>& evaluation_polynomial,
    std::vector<std::complex<double>>& lagrange_multipliers);
#endif // __cplusplus

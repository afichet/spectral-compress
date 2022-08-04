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
 * @brief Computes moments image from a spectral image
 * 
 * @param phases            Phases at which the signal was sampled.
 * @param n_phases          Number of samples.
 * @param spectral_image    Spectral image.
 * @param width             Width of the spectral image.
 * @param height            Height of the spectral image
 * @param n_moments         Number of moments to use for the moment 
 *                          representation.
 * @param moments_image     Computed moments image. Must be allocated with
 *                          `width` * `height` * (`n_moments` + 1) elements.
 */
void compute_moments_image(
    const float phases[],
    size_t n_phases,
    const float spectral_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[]);


/**
 * @brief Compress a moment image
 * 
 * @param moments_image     The moments image to be compressed.
 * @param width             Width of the moments image.
 * @param height            Height of the moments image.
 * @param n_moments         Number of moments.
 * @param compressed_moments_image Computed compressed moments image. Must be
 *                          allocated with 
 *                          `width` * `height` * (`n_moments` + 1) elements.
 */
void compress_moments_image(
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float compressed_moments_image[]);


/**
 * @brief Compress a moment image for bounded signals (reflectance)
 * 
 * @param moments_image     The moments image to be compressed.
 * @param width             Width of the moments image.
 * @param height            Height of the moments image.
 * @param n_moments         Number of moments.
 * @param compressed_moments_image Computed compressed moments image. Must be
 *                          allocated with 
 *                          `width` * `height` * (`n_moments` + 1) elements.
 */
void compress_bounded_moments_image(
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float compressed_moments_image[]);


/**
 * @brief Decompress a compressed moment image
 * 
 * @param compressed_moments_image The compressed moment image to process.
 * @param width                    Width of the compressed moment image.
 * @param height                   Height of the compressed moment image.
 * @param n_moments                Number of moments in the compressed moments
 *                                 image.
 * @param moments_image            Computed moment image. Must be allocated
 *                                 with `width` * `height` * (`n_moments` + 1)
 *                                 elements.
 */
void decompress_moments_image(
    const float compressed_moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[]);


/**
 * @brief Decompress a compressed moment image for bounded signals 
 * (reflectance)
 * 
 * @param compressed_moments_image The compressed moment image to process.
 * @param width                    Width of the compressed moment image.
 * @param height                   Height of the compressed moment image.
 * @param n_moments                Number of moments in the compressed moments
 *                                 image.
 * @param moments_image            Computed moment image. Must be allocated
 *                                 with `width` * `height` * (`n_moments` + 1)
 *                                 elements.
 */
void decompress_bounded_moments_image(
    const float compressed_moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[]);


/**
 * @brief Computes a density image matching the given moments.
 * 
 * @param phases         Phases where the density shall be computed.
 * @param n_phases       Size of phases array.
 * @param moments_image  Moments image to compute the density from.
 * @param width          Width of the moments image.
 * @param height         Height of the moments image.
 * @param n_moments      Number of moments.
 * @param density_image  Computed density image. Must be allocated with
 *                       `width` * `height` * (`n_moments` + 1) elements.
 */
void compute_density_image(
    const float phases[],
    size_t n_phases,
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float density_image[]);


/**
 * @brief Computes a bounded density image matching the given moments.
 * 
 * @param phases         Phases where the density shall be computed.
 * @param n_phases       Size of phases array.
 * @param moments_image  Moments image to compute the density from.
 * @param width          Width of the moments image.
 * @param height         Height of the moments image.
 * @param n_moments      Number of moments.
 * @param density_image  Computed density image. Must be allocated with
 *                       `width` * `height` * (`n_moments` + 1) elements.
 */
void compute_density_bounded_lagrange_image(
    const float phases[],
    size_t n_phases,
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float density_image[]);

} // extern "C"


/******************************************************************************
 * C++ interface
 *****************************************************************************/

#ifdef __cplusplus

#include <vector>

void compute_moments_image(
    const std::vector<float>& phases,
    const std::vector<float>& spectral_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& moments_image);


void compress_moments_image(
    const std::vector<float>& moments_image,
    size_t width, size_t height, 
    size_t n_moments,
    std::vector<float>& compressed_moments_image);


void decompress_moments_image(
    const std::vector<float>& compressed_moments_image,
    size_t width, size_t height, 
    size_t n_moments,
    std::vector<float>& moments_image);


void compute_density_image(
    const std::vector<float>& phases,
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& density_image);


void compute_density_bounded_lagrange_image(
    const std::vector<float>& phases,
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& density_image);

#endif // __cplusplus

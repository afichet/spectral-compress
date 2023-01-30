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

#include <spectral_compression_type.h>

#ifdef __cplusplus
    #include <cstddef>
    #include <cstdint>
    #include <cmath>
    #include <vector>
    #include <cassert>

    #include <limits>
    #include <chrono>
#else
    #include <stddef.h>
    #include <stdint.h>
    #include <math.h>
    #include <assert.h>
#endif // __cplusplus


/**
 * @brief Computes moments image from a spectral image
 *
 * @param phases            Phases at which the signal was sampled.
 * @param n_phases          Number of samples.
 * @param spectral_image    Spectral image.
 * @param n_pixels          `width * height` of the spectral image.
 * @param n_moments         Number of moments to use for the moment
 *                          representation.
 * @param moments_image     Computed moments image. Must be allocated with
 *                          `width` * `height` * `n_moments` elements.
 */
extern "C"
void compute_moments_image(
    const double phases[],
    size_t n_phases,
    const double spectral_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[]);

#ifdef __cplusplus
void compute_moments_image(
    const std::vector<double>& phases,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image);
#endif // __cplusplus


/*****************************************************************************/
/* Compression                                                               */
/*****************************************************************************/

/**
 * @brief Compress a moment image
 *
 * @param moments_image     The moments image to be compressed.
 * @param n_pixels          `width * height` of the moments image.
 * @param n_moments         Number of moments.
 * @param compressed_moments_image Computed compressed moments image. Must be
 *                          allocated with
 *                          `width` * `height` * `n_moments` elements.
 */
extern "C"
void unbounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moments_image[]);

#ifdef __cplusplus
void unbounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compressed_moments_image);
#endif // __cplusplus


/**
 * @brief Compress a moment image for bounded signals (reflectance)
 *
 * @param moments_image     The moments image to be compressed.
 * @param n_pixels          `width * height` of the moments image.
 * @param n_moments         Number of moments.
 * @param compressed_moments_image Computed compressed moments image. Must be
 *                          allocated with
 *                          `width` * `height` * `n_moments` elements.
 */
extern "C"
void bounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moments_image[]);

#ifdef __cplusplus
void bounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compressed_moments_image);
#endif // __cplusplus


/**
 * @brief Compress a moment image
 *
 * This uses zeroth moment to rescale the other moment in order to
 * use the bounded compression method
 *
 * @param moments_image     The moments image to be compressed.
 * @param n_pixels          `width * height` of the moments image.
 * @param n_moments         Number of moments.
 * @param compressed_moments_image Computed compressed moments image. Must be
 *                          allocated with
 *                          `width` * `height` * `n_moments` elements.
 */
extern "C"
void unbounded_to_bounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moments_image[]);

#ifdef __cplusplus
void unbounded_to_bounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compresed_moments_image);
#endif // __cplusplus




/*****************************************************************************/
/* Decompression                                                             */
/*****************************************************************************/

/**
 * @brief Decompress a compressed moment image
 *
 * @param compressed_moments_image The compressed moment image to process.
 * @param n_pixels                 `width * height` of the moments image.
 * @param n_moments                Number of moments in the compressed moments
 *                                 image.
 * @param moments_image            Computed moment image. Must be allocated
 *                                 with `width` * `height` * `n_moments`
 *                                 elements.
 */
extern "C"
void unbounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[]);

#ifdef __cplusplus
void unbounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image);
#endif // __cplusplus

/**
 * @brief Decompress a compressed moment image for bounded signals
 * (reflectance)
 *
 * @param compressed_moments_image The compressed moment image to process.
 * @param n_pixels                 `width * height` of the moments image.
 * @param n_moments                Number of moments in the compressed moments
 *                                 image.
 * @param moments_image            Computed moment image. Must be allocated
 *                                 with `width` * `height` * `n_moments`
 *                                 elements.
 */
extern "C"
void bounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[]);

#ifdef __cplusplus
void bounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image);
#endif // __cplusplus

/**
 * @brief Decompress a compressed moment image
 *
 * @param compressed_moments_image The compressed moment image to process.
 * @param n_pixels                 `width * height` of the moments image.
 * @param n_moments                Number of moments in the compressed moments
 *                                 image.
 * @param moments_image            Computed moment image. Must be allocated
 *                                 with `width` * `height` * `n_moments`
 *                                 elements.
 */
extern "C"
void unbounded_to_bounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[]);

#ifdef __cplusplus
void unbounded_to_bounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image);
#endif // __cplusplus

/*****************************************************************************/

/**
 * @brief Computes a density image matching the given moments.
 *
 * @param phases         Phases where the density shall be computed.
 * @param n_phases       Size of phases array.
 * @param moments_image  Moments image to compute the density from.
 * @param n_pixels       `width * height` of the moments image.
 * @param n_moments      Number of moments.
 * @param density_image  Computed density image. Must be allocated with
 *                       `width` * `height` * `n_moments` elements.
 */
extern "C"
void compute_density_image(
    const double phases[],
    size_t n_phases,
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double density_image[]);

#ifdef __cplusplus
void compute_density_image(
    const std::vector<double>& phases,
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& density_image);
#endif // __cplusplus

/**
 * @brief Computes a bounded density image matching the given moments.
 *
 * @param phases         Phases where the density shall be computed.
 * @param n_phases       Size of phases array.
 * @param moments_image  Moments image to compute the density from.
 * @param n_pixels       `width * height` of the moments image.
 * @param n_moments      Number of moments.
 * @param density_image  Computed density image. Must be allocated with
 *                       `width` * `height` * `n_moments` elements.
 */
extern "C"
void bounded_compute_density_lagrange_image(
    const double phases[],
    size_t n_phases,
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double density_image[]);

#ifdef __cplusplus
void bounded_compute_density_lagrange_image(
    const std::vector<double>& phases,
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& density_image);
#endif // __cplusplus


/**
 * @brief Rescale AC components of the image in [0, 1]
 *
 * For each AC component, rescale those to ensure all components are within
 * [0, 1]. The scaling factors are then populated in `mins[]` and `maxs[]`.
 *
 * Moments are scaled with:
 *  dest[px * n_moments + m] = (src[px * n_moments + m] - mins[m]) / (maxs[m] - mins[m])
 *
 * for all m > 0.
 *
 * @param src       Source array to rescale of size `n_pixel * n_moments`.
 * @param n_pixels  Number of the pixels in the image.
 * @param n_moments Number of moments per pixel.
 * @param dest      Destination array containig rescaled values for
 *                  moments > 0. Must be initialized with size
 *                  `n_pixel * n_moments`.
 * @param mins      Lower bound for a given moment over all pixels.
 * @param maxs      Upper bound for a given moment over all pixels.
 */
extern "C"
void normalize_moment_image_min_max(
    const double src[],
    size_t n_pixels, size_t n_moments,
    double dest[],
    double mins[],
    double maxs[]);

#ifdef __cplusplus
void normalize_moment_image_min_max(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    std::vector<double>& dest,
    std::vector<double>& mins,
    std::vector<double>& maxs);
#endif // __cplusplus


/**
 * @brief Revert a rescaed vector.
 *
 * For each AC component (moment > 0), the destination is computed as follows:
 *  dest[px * n_moments + m] = src[x * n_moments + m] * (maxs[m] - mins[m]) + mins[m]
 *
 * @param src       Source array to rescale of size `n_pixel * n_moments`.
 * @param n_pixels  Number of the pixels in the image.
 * @param n_moments Number of moments per pixel.
 * @param mins      Lower bound for a given moment over all pixels.
 * @param maxs      Upper bound for a given moment over all pixels.
 * @param dest      Destination array containig rescaled values for
 *                  moments > 0. Must be initialized with size
 *                  `n_pixel * n_moments`.
 */
extern "C"
void denormalize_moment_image_min_max(
    const double src[],
    size_t n_pixels, size_t n_moments,
    const double mins[],
    const double maxs[],
    double dest[]);

#ifdef __cplusplus
void denormalize_moment_image_min_max(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    std::vector<double>& dest);
#endif // __cplusplus


/**
 * @brief Rescale AC components of the image in [0, 1]
 *
 * For each AC component, rescale those to ensure all components are within
 * [0, 1]. The scaling factors are then populated in `mins[]` and `maxs[]`.
 *
 * Moments are scaled with:
 *  dest[px * n_moments + m] = (src[px * n_moments + m] - mins[m]) / (maxs[m] - mins[m])
 *
 * for all m > 0.
 *
 * @param src       Source array to rescale of size `n_pixel * n_moments`.
 * @param n_pixels  Number of the pixels in the image.
 * @param n_moments Number of moments per pixel.
 * @param dest      Destination array containig rescaled values for
 *                  moments > 0. Must be initialized with size
 *                  `n_pixel * n_moments`.
 * @param means     Mean for a given moment over all pixels.
 * @param stddevs   Standard deviation for a given moment over all pixels.
 */
extern "C"
void normalize_moment_image_stddev(
    const double src[],
    size_t n_pixels, size_t n_moments,
    double dest[],
    double means[],
    double stddevs[]);

#ifdef __cplusplus
void normalize_moment_image_stddev(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    std::vector<double>& dest,
    std::vector<double>& means,
    std::vector<double>& stddevs);
#endif // __cplusplus


/**
 * @brief Revert a rescaed vector.
 *
 * For each AC component (moment > 0), the destination is computed as follows:
 *  dest[px * n_moments + m] = src[x * n_moments + m] * (maxs[m] - mins[m]) + mins[m]
 *
 * @param src       Source array to rescale of size `n_pixel * n_moments`.
 * @param n_pixels  Number of the pixels in the image.
 * @param n_moments Number of moments per pixel.
 * @param means     Mean for a given moment over all pixels.
 * @param stddevs   Standard deviation for a given moment over all pixels.
 * @param dest      Destination array containig rescaled values for
 *                  moments > 0. Must be initialized with size
 *                  `n_pixel * n_moments`.
 */
extern "C"
void denormalize_moment_image_stddev(
    const double src[],
    size_t n_pixels, size_t n_moments,
    const double means[],
    const double stddevs[],
    double dest[]);

#ifdef __cplusplus
void denormalize_moment_image_stddev(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    const std::vector<double>& means,
    const std::vector<double>& stddevs,
    std::vector<double>& dest);
#endif // __cplusplus


/*****************************************************************************/
/* Full pipeline                                                             */
/*****************************************************************************/

// Compress (spectral_image -> compressed_moments)

#ifdef __cplusplus

void linear_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs,
    bool normalize_image);


void linavg_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs,
    bool normalize_image);


// All in one function
void compress_spectral_image(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double> spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    std::vector<double>& compressed_moments,
    std::vector<double>& mins, std::vector<double>& maxs,
    bool normalize_image,
    double& timing);


// Decompress (compressed_moments -> spectral_image)

void linear_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image);


void linavg_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image);


// All in one function
void decompress_spectral_image(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image,
    double& timing);

#endif // __cplusplus

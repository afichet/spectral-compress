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

#include <vector>
#include <cassert>
#include <cstdint>
#include <cstddef>

enum CompressionCurveType {
    COMPRESSION_FLAT,
    COMPRESSION_DYNAMIC,
    COMPRESSION_DETERMINISTIC
};

/*****************************************************************************/
/* Utility functions                                                         */
/*****************************************************************************/

void compress_decompress_framebuffer(
    const std::vector<float>& framebuffer_in,
    std::vector<float>& framebuffer_out,
    uint32_t width, uint32_t height,
    std::pair<int, int> bits_per_sample,
    uint32_t subsampling_factor,
    float frame_distance,
    int effort);


void compress_decompress_single_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    uint32_t bits_per_sample,
    uint32_t subsampling_factor,
    size_t i, float frame_distance,
    int effort);


void compress_decompress_image(
    const std::vector<double>& input_framebuffers,
    std::vector<double>& output_framebuffers,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


/*****************************************************************************/
/* Create compression curves                                                 */
/*****************************************************************************/

void compute_compression_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    CompressionCurveType compression_curve_type,
    std::vector<float>& compression_curve,
    float effort,
    double& timing);


double linear_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double linavg_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double unbounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double unbounded_to_bounded_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double upperbound_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


double twobounds_compute_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    bool normalize_moments,
    const std::vector<uint32_t>& subsampling_factor_curve,
    float compression_dc,
    float compression_ac1,
    std::vector<float>& compression_curve,
    int effort = 7);


/*****************************************************************************/
/* Error for a compression curve                                             */
/*****************************************************************************/

double linear_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double linavg_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double unbounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double bounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double unbounded_to_bounded_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double upperbound_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


double twobounds_error_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_min,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& subsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);

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
#include <Util.h>

#include <vector>
#include <cstdint>
#include <cstddef>

struct stats_data {
    double rmse_error;
    double rrmse_error;
    double avg_error;
    double max_error;
};

template<typename T>
stats_data compute_stats(
    const std::vector<T>& reference,
    const std::vector<T>& comparison,
    uint32_t width, uint32_t height,
    size_t n_wavelengths)
{
        stats_data ret;

    ret.rmse_error = Util::rmse_images(
        reference,
        comparison,
        width * height,
        n_wavelengths
    );

    ret.rrmse_error = Util::rrmse_images(
        reference,
        comparison,
        width * height,
        n_wavelengths
    );

    ret.avg_error = Util::avg_err_images(
        reference,
        comparison,
        width * height,
        n_wavelengths
    );

    ret.max_error = Util::max_error_images(
        reference,
        comparison,
        width * height,
        n_wavelengths
    );

    return ret;
}

/*****************************************************************************/
/* Stats for quantization curves                                             */
/*****************************************************************************/

// More or less the same thing as in `common/curve_quantization.h` but
// with more statistics

stats_data stats_for_quantization_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max);


stats_data linear_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve);


stats_data linavg_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve);


stats_data unbounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve);


stats_data bounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve);


stats_data unbounded_to_bounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve);


// TODO: template?
stats_data upperbound_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve);


// TODO: template?
stats_data twobounds_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_min,
    double global_max,
    const std::vector<std::pair<int, int>>& quantization_curve);


/*****************************************************************************/
/* Stats for compression curves                                              */
/*****************************************************************************/

// More or less the same thing as in `common/curve_compression.h` but
// with more statistics

stats_data stats_for_compression_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max,
    int effort);


stats_data linear_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data linavg_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data unbounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data bounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data unbounded_to_bounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::pair<int, int>>& quantization_curve,
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data upperbound_stats_for_compression_curve(
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
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);


stats_data twobounds_stats_for_compression_curve(
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
    const std::vector<uint32_t>& downsampling_factor_curve,
    const std::vector<float>& compression_curve,
    int effort);

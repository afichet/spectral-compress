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

#include <vector>
#include <cstddef>

#include "Util.h"


double linear_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs);


double unbounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs);


double bounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs);


double unbounded_to_bounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs);


template<typename T>
double upperbound_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<T>& relative_scales,
    double global_max)
{
    // double err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<double> reconst_spectral_image;

    upperbound_decompress_spectral_image(
        wavelengths, norm_moments,
        mins, maxs,
        relative_scales,
        global_max,
        n_px, n_moments,
        reconst_spectral_image
    );

    return Util::rmse_images(
        spectral_image,
        reconst_spectral_image,
        n_px, n_wl
    );
}


template<typename T>
double twobounds_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<T>& relative_scales,
    double global_min,
    double global_max)
{
    // double err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<double> reconst_spectral_image;

    twobounds_decompress_spectral_image(
        wavelengths, norm_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        n_px, n_moments,
        reconst_spectral_image
    );

    return Util::rmse_images(
        spectral_image,
        reconst_spectral_image,
        n_px, n_wl
    );
}

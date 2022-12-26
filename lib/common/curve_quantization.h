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

#include <vector>
#include <cstdint>
#include <cstddef>


/*****************************************************************************/
/* Utility functions                                                         */
/*****************************************************************************/

void quantize_dequantize_single_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    size_t n_pixels, size_t n_moments,
    size_t i, size_t n_bits);


void quantize_dequantize_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    size_t n_pixels, size_t n_moments,
    const std::vector<int>& quantization_curve);


/*****************************************************************************/
/* Create quantization curves                                                */
/*****************************************************************************/

double linear_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


double unbounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


double bounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


double unbounded_to_bounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


double upperbound_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


double twobounds_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve);


/*****************************************************************************/
/* Error for a quantization curve                                            */
/*****************************************************************************/

double linear_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve);


double unbounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve);


double bounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve);


double unbounded_to_bounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve);


// TODO: template?
double upperbound_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_max,
    const std::vector<int>& quantization_curve);


// TODO: template?
double twobounds_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_min,
    double global_max,
    const std::vector<int>& quantization_curve);

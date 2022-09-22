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

#include "quantization.h"
#include "moments.h"
#include "moments_image.h"

#include <cmath>
#include <cstring>

// #define USE_BOUNDED

inline float quantize_dequantize(float src, size_t n_bits)
{
    const uint64_t v = (1 << n_bits) - 1;

    return std::round(src * (float)v) / (float)v; 
}


void quantize_dequantize_single_image(
    const std::vector<float>& src,
    size_t n_px, size_t n_moments,
    std::vector<float>& dest,
    size_t i, size_t n_bits)
{
    dest.resize(src.size());

    std::memcpy(dest.data(), src.data(), sizeof(float) * src.size());

    for (size_t px = 0; px < n_px; px++) {
        dest[px * n_moments + i] = quantize_dequantize(src[px * n_moments + i], n_bits);
    }
}


float average_err(
    const std::vector<float>& wavelengths,
    const std::vector<float>& ref,
    size_t n_px, size_t n_moments,
    const std::vector<float>& norm_noments,
    const std::vector<float>& mins,
    const std::vector<float>& maxs)
{
    float err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<float> phases;

    std::vector<float> compressed_moments;
    std::vector<float> moments;
    std::vector<float> reconst_signal;

    wavelengths_to_phases(wavelengths, phases);
    denormalize_moment_image(norm_noments, n_px, n_moments, mins, maxs, compressed_moments);
#if defined(USE_BOUNDED)
    decompress_bounded_moments_image(compressed_moments, n_px, 1, n_moments, moments);
#else
    decompress_moments_image(compressed_moments, n_px, 1, n_moments, moments);
#endif // USE_BOUNDED
    compute_density_image(phases, moments, n_px, 1, n_moments, reconst_signal);

    for (size_t p = 0; p < n_px; p++) {
        for (size_t i = 0; i < n_wl; i++) {
            const float q = reconst_signal[p * n_wl + i] - ref[p * n_wl + i];
            err += q * q;
        }
    }

    return err;
}


void compute_quantization_curve(
    const std::vector<float>& wavelengths,
    const std::vector<float>& ref,
    size_t n_px, size_t n_moments,
    const std::vector<float>& norm_moments,
    const std::vector<float>& mins,
    const std::vector<float>& maxs,
    int n_bits_start,
    std::vector<int>& quantization_curve)
{
    quantization_curve.resize(n_moments);

    quantization_curve[0] = 16; // TODO: remove
    quantization_curve[1] = n_bits_start;

    std::vector<float> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        n_px, n_moments,
        quantized_moments,
        1, quantization_curve[1]
    );

    const float base_err = average_err(
        wavelengths,
        ref,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );


    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m]; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                n_px, n_moments,
                quantized_moments,
                m, n_bits
            );

            float curr_err = average_err(
                wavelengths,
                ref,
                n_px, n_moments,
                quantized_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            quantization_curve[m] = n_bits;
        }
    }
}
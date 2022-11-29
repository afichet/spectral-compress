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
#include <cassert>


inline double quantize_dequantize(double src, size_t n_bits)
{
    const uint64_t v = (1 << n_bits) - 1;

    return std::round(src * (double)v) / (double)v;
}


void quantize_dequantize_single_image(
    const std::vector<double>& src,
    size_t n_px, size_t n_moments,
    std::vector<double>& dest,
    size_t i, size_t n_bits)
{
    dest.resize(src.size());

    std::memcpy(dest.data(), src.data(), sizeof(double) * src.size());

    for (size_t px = 0; px < n_px; px++) {
        dest[px * n_moments + i] = quantize_dequantize(src[px * n_moments + i], n_bits);
    }
}


double unbounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs)
{
    double err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<double> reconst_spectral_image;

    unbounded_decompress_spectral_image(
        wavelengths, norm_moments,
        mins, maxs,
        n_px, n_moments,
        reconst_spectral_image
    );

    for (size_t p = 0; p < n_px; p++) {
        double px_err = 0;

        if (norm_moments[p * n_wl] > 0) {
            for (size_t i = 0; i < n_wl; i++) {
                const double q = reconst_spectral_image[p * n_wl + i] - spectral_image[p * n_wl + i];
                px_err += q * q;
            }

            px_err = std::sqrt(px_err) / norm_moments[p * n_wl];
        }

        err += px_err;
    }

    return err / (double)n_px;
}


double bounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs)
{
    double err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<double> reconst_spectral_image;

    bounded_decompress_spectral_image(
        wavelengths, norm_moments,
        mins, maxs,
        n_px, n_moments,
        reconst_spectral_image
    );

    for (size_t p = 0; p < n_px; p++) {
        double px_err = 0;

        if (norm_moments[p * n_wl] > 0) {
            for (size_t i = 0; i < n_wl; i++) {
                const double q = reconst_spectral_image[p * n_wl + i] - spectral_image[p * n_wl + i];
                px_err += q * q;
            }

            px_err = std::sqrt(px_err) / norm_moments[p * n_wl];
        }

        err += px_err;
    }

    return err / (double)n_px;
}


double unbounded_to_bounded_average_err(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& norm_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs)
{
    double err = 0;
    const size_t n_wl = wavelengths.size();

    std::vector<double> reconst_spectral_image;

    unbounded_to_bounded_decompress_spectral_image(
        wavelengths, norm_moments,
        mins, maxs,
        n_px, n_moments,
        reconst_spectral_image
    );

    for (size_t p = 0; p < n_px; p++) {
        double px_err = 0;

        if (norm_moments[p * n_wl] > 0) {
            for (size_t i = 0; i < n_wl; i++) {
                const double q = reconst_spectral_image[p * n_wl + i] - spectral_image[p * n_wl + i];
                px_err += q * q;
            }

            px_err = std::sqrt(px_err) / norm_moments[p * n_wl];
        }

        err += px_err;
    }

    return err / (double)n_px;
}


/*****************************************************************************/
/* Create quantization curves                                                */
/*****************************************************************************/

double unbounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_start,
    std::vector<int>& quantization_curve,
    int n_bits_0)
{
    std::vector<double> norm_moments;
    std::vector<double> mins, maxs;

    unbounded_compress_spectral_image(
        wavelengths, spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_0;
    quantization_curve[1] = n_bits_start;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        n_px, n_moments,
        quantized_moments,
        1, quantization_curve[1]
    );

    const double base_err = unbounded_average_err(
        wavelengths,
        spectral_image,
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

            double curr_err = unbounded_average_err(
                wavelengths,
                spectral_image,
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

    return unbounded_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs,
        quantization_curve
    );
}


double bounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_start,
    std::vector<int>& quantization_curve,
    int n_bits_0)
{
    std::vector<double> norm_moments;
    std::vector<double> mins, maxs;

    bounded_compress_spectral_image(
        wavelengths, spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_0;
    quantization_curve[1] = n_bits_start;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        n_px,
        n_moments,
        quantized_moments,
        1, quantization_curve[1]
    );

    const double base_err = bounded_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );


    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m]; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                n_px,
                n_moments,
                quantized_moments,
                m, n_bits
            );

            double curr_err = bounded_average_err(
                wavelengths,
                spectral_image,
                n_px,
                n_moments,
                quantized_moments,
                mins, maxs
            );

            if (curr_err >= base_err) {
                break;
            }

            quantization_curve[m] = n_bits;
        }
    }

    return bounded_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px,
        n_moments,
        norm_moments,
        mins, maxs,
        quantization_curve
    );
}


double unbounded_to_bounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_start,
    std::vector<int>& quantization_curve,
    int n_bits_0)
{
    std::vector<double> norm_moments;
    std::vector<double> mins, maxs;

    unbounded_to_bounded_compress_spectral_image(
        wavelengths, spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_0;
    quantization_curve[1] = n_bits_start;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        n_px, n_moments,
        quantized_moments,
        1, quantization_curve[1]
    );

    const double base_err = unbounded_to_bounded_average_err(
        wavelengths,
        spectral_image,
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

            double curr_err = unbounded_to_bounded_average_err(
                wavelengths,
                spectral_image,
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


    return unbounded_to_bounded_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs,
        quantization_curve
    );
}


/*****************************************************************************/
/* Error for a quantization curve                                            */
/*****************************************************************************/

double unbounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < n_px; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    return unbounded_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}


double bounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < n_px; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    return bounded_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}


double unbounded_to_bounded_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments(normalized_moments.size());

    for (size_t px = 0; px < n_px; px++) {
        // We ignore moment 0
        quantized_moments[px * n_moments + 0] = normalized_moments[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            quantized_moments[px * n_moments + m] =
                quantize_dequantize(
                    normalized_moments[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }

    return unbounded_to_bounded_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}

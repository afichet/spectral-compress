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

#include "curve_quantization.h"
#include "moments.h"
#include "moments_image.h"
#include "moments_error.h"

#include <cmath>
#include <cstring>
#include <cassert>
#include <chrono>


void quantize_dequantize_single_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    size_t n_pixels, size_t n_moments,
    size_t i, size_t n_bits)
{
    assert(input_image.size() == n_pixels * n_moments);
    assert(i < n_moments);

    output_image.resize(input_image.size());

    std::memcpy(output_image.data(), input_image.data(), sizeof(double) * input_image.size());

    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        output_image[px * n_moments + i] = Util::quantize_dequantize(input_image[px * n_moments + i], n_bits);
    }
}


void quantize_dequantize_image(
    const std::vector<double>& input_image,
    std::vector<double>& output_image,
    size_t n_pixels, size_t n_moments,
    const std::vector<int>& quantization_curve)
{
    assert(input_image.size() == n_pixels * n_moments);
    assert(quantization_curve.size() >= n_moments);

    output_image.resize(input_image.size());

    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // We ignore moment 0
        output_image[px * n_moments + 0] = input_image[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            output_image[px * n_moments + m] =
                Util::quantize_dequantize(
                    input_image[px * n_moments + m],
                    quantization_curve[m]
                );
        }
    }
}

/*****************************************************************************/
/* Create quantization curves                                                */
/*****************************************************************************/

void compute_quantization_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    bool uses_constant_quantization,
    std::vector<int>& quantization_curve,
    double& timing)
{
    if (uses_constant_quantization) {
        quantization_curve.resize(n_moments);
        quantization_curve[0] = n_bits_dc;

        for (size_t i = 1; i < n_moments; i++) {
            quantization_curve[i] = n_bits_ac1;
        }
        timing = 0;
    } else {
        auto clock_start = std::chrono::steady_clock::now();

        switch (method) {
            case LINEAR:
                linear_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case LINAVG:
                linavg_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case BOUNDED:
                bounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UNBOUNDED:
                unbounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UNBOUNDED_TO_BOUNDED:
                unbounded_to_bounded_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case UPPERBOUND:
                upperbound_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
            case TWOBOUNDS:
                twobounds_compute_quantization_curve(
                    wavelengths, spectral_image,
                    width * height, n_moments,
                    n_bits_dc, n_bits_ac1,
                    quantization_curve
                );
                break;
        }

        auto clock_end = std::chrono::steady_clock::now();
        timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();
    }

    assert(quantization_curve.size() == n_moments);
}


double linear_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
{
    std::vector<double> norm_moments;
    std::vector<double> mins, maxs;

    linear_compress_spectral_image(
        wavelengths, spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        quantized_moments,
        n_px, n_moments,
        1, quantization_curve[1]
    );

    const double base_err = linear_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = linear_average_err(
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

    return linear_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs,
        quantization_curve
    );
}


double linavg_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
{
    std::vector<double> norm_moments;
    std::vector<double> mins, maxs;

    linavg_compress_spectral_image(
        wavelengths, spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        quantized_moments,
        n_px, n_moments,
        1, quantization_curve[1]
    );

    const double base_err = linavg_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );

    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = linavg_average_err(
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

    return linavg_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        norm_moments,
        mins, maxs,
        quantization_curve
    );
}


double unbounded_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
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

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        quantized_moments,
        n_px, n_moments,
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

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = unbounded_average_err(
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
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
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

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        quantized_moments,
        n_px, n_moments,
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

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = bounded_average_err(
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
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
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

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        norm_moments,
        quantized_moments,
        n_px, n_moments,
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

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                norm_moments,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = unbounded_to_bounded_average_err(
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

// TODO: use template?
double upperbound_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
{
    std::vector<double> normalized_moments_image;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_max;

    upperbound_compress_spectral_image(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        normalized_moments_image,
        mins, maxs,
        relative_scales,
        global_max
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        normalized_moments_image,
        quantized_moments,
        n_px, n_moments,
        1, quantization_curve[1]
    );

    const double base_err = upperbound_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs,
        relative_scales,
        global_max
    );

    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                normalized_moments_image,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = upperbound_average_err(
                wavelengths,
                spectral_image,
                n_px, n_moments,
                quantized_moments,
                mins, maxs,
                relative_scales,
                global_max
            );

            if (curr_err >= base_err) {
                break;
            }

            quantization_curve[m] = n_bits;
        }
    }

    return upperbound_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        normalized_moments_image,
        mins, maxs,
        relative_scales,
        global_max,
        quantization_curve
    );
}


// TODO: use template?
double twobounds_compute_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    int n_bits_dc,
    int n_bits_ac1,
    std::vector<int>& quantization_curve)
{
    std::vector<double> normalized_moments_image;
    std::vector<double> mins, maxs;
    std::vector<uint8_t> relative_scales;
    double global_min, global_max;

    twobounds_compress_spectral_image(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        normalized_moments_image,
        mins, maxs,
        relative_scales,
        global_min,
        global_max
    );

    quantization_curve.resize(n_moments);

    quantization_curve[0] = n_bits_dc;
    quantization_curve[1] = n_bits_ac1;

    std::vector<double> quantized_moments;

    quantize_dequantize_single_image(
        normalized_moments_image,
        quantized_moments,
        n_px, n_moments,
        1, quantization_curve[1]
    );

    const double base_err = twobounds_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max
    );

    for (size_t m = 2; m < n_moments; m++) {
        quantization_curve[m] = quantization_curve[m - 1];

        for (size_t n_bits = quantization_curve[m] - 1; n_bits > 0; n_bits--) {
            quantize_dequantize_single_image(
                normalized_moments_image,
                quantized_moments,
                n_px, n_moments,
                m, n_bits
            );

            const double curr_err = twobounds_average_err(
                wavelengths,
                spectral_image,
                n_px, n_moments,
                quantized_moments,
                mins, maxs,
                relative_scales,
                global_min,
                global_max
            );

            if (curr_err >= base_err) {
                break;
            }

            quantization_curve[m] = n_bits;
        }
    }

    return twobounds_error_for_quantization_curve(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        normalized_moments_image,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        quantization_curve
    );
}


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
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

    return linear_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}


double linavg_error_for_quantization_curve(
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
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

    return linavg_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}


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
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

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
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

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
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

    return unbounded_to_bounded_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs
    );
}


double upperbound_error_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double global_max,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

    return upperbound_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs,
        relative_scales,
        global_max
    );
}


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
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == n_px * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        n_px, n_moments,
        quantization_curve
    );

    return twobounds_average_err(
        wavelengths,
        spectral_image,
        n_px, n_moments,
        quantized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max
    );
}

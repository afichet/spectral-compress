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


#include "stats.h"

#include <cassert>

#include <curve_quantization.h>
#include <curve_compression.h>
#include <moments_image.h>

#include <Util.h>


/*****************************************************************************/
/* Stats for quantization curves                                             */
/*****************************************************************************/

stats_data linear_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    linear_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        width * height,
        n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data unbounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    unbounded_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        width * height,
        n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data bounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    bounded_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        width * height,
        n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data unbounded_to_bounded_stats_for_quantization_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    unbounded_to_bounded_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        width * height,
        n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


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
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height,
        n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    upperbound_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        relative_scales,
        global_max,
        width * height,
        n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


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
    const std::vector<int>& quantization_curve)
{
    assert(normalized_moments.size() == width * height * n_moments);

    // Quantize & dequantize each moment of the image
    std::vector<double> quantized_moments;

    quantize_dequantize_image(
        normalized_moments,
        quantized_moments,
        width * height, n_moments,
        quantization_curve
    );

    std::vector<double> reconst_spectral_image;

    twobounds_decompress_spectral_image(
        wavelengths, quantized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        width * height, n_moments,
        reconst_spectral_image
    );

    return compute_stats(
        spectral_image,
        reconst_spectral_image,
        width, height,
        wavelengths.size()
    );
}


/*****************************************************************************/
/* Stats for compression curves                                              */
/*****************************************************************************/

stats_data linear_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    linear_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data unbounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    unbounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data bounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    bounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}


stats_data unbounded_to_bounded_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    unbounded_to_bounded_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}


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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    upperbound_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        relative_scales,
        global_max,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}


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
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    compress_decompress_image(
        normalized_moments,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    twobounds_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins, maxs,
        relative_scales,
        global_min,
        global_max,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    return compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );
}

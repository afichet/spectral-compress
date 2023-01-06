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

#include <curve_quantization.h>
#include <curve_compression.h>
#include <moments_image.h>

#include <cassert>

/*****************************************************************************/
/* Stats for quantization curves                                             */
/*****************************************************************************/

stats_data stats_for_quantization_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max)
{
    switch(method) {
        case LINEAR:
            return linear_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case BOUNDED:
            return bounded_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UNBOUNDED:
            return unbounded_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UNBOUNDED_TO_BOUNDED:
            return unbounded_to_bounded_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                quantization_curve
            );
        case UPPERBOUND:
            return upperbound_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_max,
                quantization_curve
            );
        case TWOBOUNDS:
            return twobounds_stats_for_quantization_curve(
                wavelengths, spectral_image,
                width, height,
                n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_min, global_max,
                quantization_curve
            );
    }

    assert(0);
    return stats_data();
}


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

stats_data stats_for_compression_curve(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve,
    const std::vector<double>& compressed_moments,
    const std::vector<double>& mins, const std::vector<double>& maxs,
    const std::vector<uint8_t>& relative_scales,
    double& global_min, double& global_max,
    int effort)
{
    switch(method) {
        case LINEAR:
            return linear_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve,
                effort
            );
        case BOUNDED:
            return bounded_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve,
                effort
            );
        case UNBOUNDED:
            return unbounded_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve,
                effort
            );
        case UNBOUNDED_TO_BOUNDED:
            return unbounded_to_bounded_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                quantization_curve,
                compression_curve,
                effort
            );
        case UPPERBOUND:
            return upperbound_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_max,
                quantization_curve,
                compression_curve,
                effort
            );
        case TWOBOUNDS:
            return twobounds_stats_for_compression_curve(
                wavelengths, spectral_image,
                width, height, n_moments,
                compressed_moments, mins, maxs,
                relative_scales, global_min, global_max,
                quantization_curve,
                compression_curve,
                effort
            );
    }

    assert(0);
    return stats_data();
}


stats_data linear_stats_for_compression_curve(
    const std::vector<double>& wavelengths,
    const std::vector<double>& ref_spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<int>& quantization_curve,
    const std::vector<float>& compression_curve,
    int effort)
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
        compression_curve,
        effort
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
    const std::vector<float>& compression_curve,
    int effort)
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
        compression_curve,
        effort
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
    const std::vector<float>& compression_curve,
    int effort)
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
        compression_curve,
        effort
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
    const std::vector<float>& compression_curve,
    int effort)
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
        compression_curve,
        effort
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
    const std::vector<float>& compression_curve,
    int effort)
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
        compression_curve,
        effort
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
    const std::vector<float>& compression_curve,
    int effort)
{
    assert(ref_spectral_image.size() == width * height * wavelengths.size());
    assert(normalized_moments.size() == width * height * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
    assert(relative_scales.size() == width * height);
    assert(quantization_curve.size() == n_moments);
    assert(compression_curve.size() == n_moments);

    std::vector<double> compressed_decompressed_normalized_moments;

    // compress_decompress_image(
    //     normalized_moments,
    //     compressed_decompressed_normalized_moments,
    //     width, height,
    //     n_moments,
    //     quantization_curve,
    //     compression_curve,
    //     effort
    // );

    // std::vector<double> decompressed_spectral_image;

    // twobounds_decompress_spectral_image(
    //     wavelengths,
    //     compressed_decompressed_normalized_moments,
    //     mins, maxs,
    //     relative_scales,
    //     global_min,
    //     global_max,
    //     width * height,
    //     n_moments,
    //     decompressed_spectral_image
    // );

    // return compute_stats(
    //     ref_spectral_image,
    //     decompressed_spectral_image,
    //     width, height,
    //     wavelengths.size()
    // );

    // TODO remove
    std::vector<float> normalized_moments_f;
    std::vector<double> normalized_moments_d;

    Util::cast_vector(normalized_moments, normalized_moments_f);
    Util::cast_vector(normalized_moments_f, normalized_moments_d);

    compress_decompress_image(
        normalized_moments_d,
        compressed_decompressed_normalized_moments,
        width, height,
        n_moments,
        quantization_curve,
        compression_curve,
        effort
    );

    // Unpack the moments
    std::vector<double> decompressed_spectral_image;

    // TODO: remove this mess
    std::vector<float> mins_f, maxs_f;
    std::vector<double> mins_d, maxs_d;

    Util::cast_vector(mins, mins_f);
    Util::cast_vector(maxs, maxs_f);

    Util::cast_vector(mins_f, mins_d);
    Util::cast_vector(maxs_f, maxs_d);
    float g_min = global_min;
    float g_max = global_max;

    twobounds_decompress_spectral_image(
        wavelengths,
        compressed_decompressed_normalized_moments,
        mins_d, maxs_d,
        relative_scales,
        (double)g_min,
        (double)g_max,
        width * height,
        n_moments,
        decompressed_spectral_image
    );

    stats_data st = compute_stats(
        ref_spectral_image,
        decompressed_spectral_image,
        width, height,
        wavelengths.size()
    );

    return st;
}

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

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <array>

#include <moments.h>
#include <moments_image.h>
#include <quantization.h>
#include <Util.h>
#include <SpectrumConverter.h>

#include "macbeth_data.h"

#define CURVES
// #define PARETO


// TODO: sorry that's ugly, I know...
SpectrumConverter sc;


int main(int argc, char* argv[])
{
    (void)argc; (void)argv;
    
    const size_t n_wl = 36;
    const size_t n_px = 24;
    const size_t n_moments = n_wl;

    std::vector<float> wavelengths(n_wl);
    for (size_t i = 0; i < n_wl; i++) {
        wavelengths[i] = macbeth_wavelengths[i];
    }

    std::vector<float> macbeth(n_wl * n_px);
    for (size_t wl = 0; wl < n_wl; wl++) {
        for (size_t px = 0; px < n_px; px++) {
            macbeth[px * n_wl + wl] = macbeth_patches[px][wl];
        }
    }

    std::vector<float> phases(n_wl);

    wavelengths_to_phases(macbeth_wavelengths, phases.size(), phases.data());

    std::vector<float> moments(n_wl * n_px);
    std::vector<float> compressed_moments(n_wl * n_px);
    std::vector<float> normalized_moments(n_wl * n_px);
    std::vector<float> mins, maxs;

    compute_moments_image(
        phases, macbeth,
        n_px, 1, n_moments,
        moments
    );

#if defined(USE_BOUNDED)
    bounded_compress_moments_image(
        moments,
        n_px, 1, n_moments,
        compressed_moments
    );
#else
    unbounded_compress_moments_image(
        moments,
        n_px, 1, n_moments,
        compressed_moments
    );
#endif

    normalize_moment_image(
        compressed_moments,
        n_px, n_moments,
        normalized_moments,
        mins, maxs
    );

    // ------------------------------------------------------------------------
    // Eval quantization
    // ------------------------------------------------------------------------

    std::vector<float> quantized_moments;

#if defined(MOMENT)
    // ------------------------------------------------------------------------
    // Error for quantizing a given moment
    // ------------------------------------------------------------------------

    for (int m = 1; m < n_moments; m++) {
        quantize_dequantize_single_image(normalized_moments, n_px, n_moments, quantized_moments, m, 8);
        const float err = average_err(wavelengths, macbeth, n_px, n_moments, quantized_moments, mins, maxs);
        std::cout << m << " " << err << std::endl;
    }

#elif defined(CURVES)
    // ------------------------------------------------------------------------
    // Determining quantization curves
    // ------------------------------------------------------------------------

    for (size_t n_bits = 3; n_bits <= 12; n_bits++) {
        std::vector<int> quantization_curve;

        compute_quantization_curve(
            wavelengths,
            macbeth,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            n_bits,
            quantization_curve
        );

        for (size_t m = 1; m < n_moments; m++) {
            std::cout << m << " " << quantization_curve[m] << std::endl;
        }

        std::cout << std::endl << std::endl;
    }

#elif defined(PARETO)
    // ------------------------------------------------------------------------
    // Pareto
    // ------------------------------------------------------------------------

    // Optimized curves
    std::vector<std::vector<int>> optimized_quantization_curves;

    for (size_t n_bits = 3; n_bits <= 12; n_bits++) {
        // Get the "optimal quantization curve"
        std::vector<int> quantization_curve;

        compute_quantization_curve(
            wavelengths,
            macbeth,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            n_bits,
            quantization_curve
        );

        size_t total_bits = 32; // starting at 1, first moment is 32 bits float

        for (size_t m = 1; m < n_moments; m++) {
            total_bits += quantization_curve[m];    
        }

        float error = error_for_quantization_curve(
            wavelengths, macbeth,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            quantization_curve
        );

        std::cout << total_bits << " " << error << std::endl;

        optimized_quantization_curves.push_back(quantization_curve);
    }

    std::cout << std::endl << std::endl;

    // Mutated curves
    for (size_t q_idx = 0; q_idx < optimized_quantization_curves.size(); q_idx++) {
        // Mutate the curve
        const int n_mutations = 1000;

        for (size_t i = 0; i < n_mutations; i++) {
            size_t total_bits = 32; // starting at 1, first moment is 32 bits float

            std::vector<int> mutated_curve = optimized_quantization_curves[q_idx];

            for (size_t m = 1; m < n_moments; m++) {
                const int mutation = std::round(16.f * ((float)std::rand() / (float)RAND_MAX) - 8.f);
                mutated_curve[m] += std::max(mutated_curve[m] + mutation, 1);

                total_bits += mutated_curve[m];
            }

            const float error = error_for_quantization_curve(
                wavelengths, macbeth,
                n_px, n_moments,
                normalized_moments,
                mins, maxs,
                mutated_curve
            );

            std::cout << total_bits << " " << error << std::endl;
        }
    }

    std::cout << std::endl << std::endl;

    // Randomized curves
    for (size_t i = 0; i < 1000; i++) {
        for (size_t n_bits = 2; n_bits <= 12; n_bits++) {
            // Generate a random quantization curve
            std::vector<int> quantization_curve(n_moments);
            
            size_t total_bits = 32; // starting at 1, first moment is 32 bits float

            // WARN: we start at 1
            for (size_t m = 1; m < n_moments; m++) {
                quantization_curve[m] = std::round((n_bits - 1) * (float)std::rand() / (float)RAND_MAX) + 1;

                total_bits += quantization_curve[m];
            }

            const float error = error_for_quantization_curve(
                wavelengths, macbeth,
                n_px, n_moments,
                normalized_moments,
                mins, maxs,
                quantization_curve
            );

            std::cout << total_bits << " " << error << std::endl;
        }
    }
#endif

    return 0;
}
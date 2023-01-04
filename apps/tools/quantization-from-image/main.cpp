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

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <array>
#include <fstream>

#include <moments.h>
#include <moments_image.h>
#include <curve_quantization.h>
#include <Util.h>
#include <SpectrumConverter.h>
#include <EXRSpectralImage.h>

#include "macbeth_data.h"

#define CURVES
// #define PARETO


// TODO: sorry that's ugly, I know...
SpectrumConverter sc;

// ----------------------------------------------------------------------------
// Optimal quantization curve generation
// ----------------------------------------------------------------------------

void bounded_generate_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    std::vector<std::vector<int>>& optimized_quantization_curves)
{
    optimized_quantization_curves.clear();

    for (size_t n_bits = 8; n_bits <= 12; n_bits++) {
        optimized_quantization_curves.resize(optimized_quantization_curves.size() + 1);

        std::vector<int>& quantization_curve = optimized_quantization_curves.back();

        bounded_compute_quantization_curve(
            wavelengths,
            spectral_image,
            n_px, n_moments,
            32,
            n_bits,
            quantization_curve
        );
    }
}


void unbounded_generate_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    std::vector<std::vector<int>>& optimized_quantization_curves)
{
    optimized_quantization_curves.clear();

    for (size_t n_bits = 8; n_bits <= 12; n_bits++) {
        optimized_quantization_curves.resize(optimized_quantization_curves.size() + 1);

        std::vector<int>& quantization_curve = optimized_quantization_curves.back();

        unbounded_compute_quantization_curve(
            wavelengths,
            spectral_image,
            n_px, n_moments,
            32,
            n_bits,
            quantization_curve
        );
    }
}


void unbounded_to_bounded_generate_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    std::vector<std::vector<int>>& optimized_quantization_curves)
{
    optimized_quantization_curves.clear();

    for (size_t n_bits = 8; n_bits <= 12; n_bits++) {
        optimized_quantization_curves.resize(optimized_quantization_curves.size() + 1);

        std::vector<int>& quantization_curve = optimized_quantization_curves.back();

        unbounded_to_bounded_compute_quantization_curve(
            wavelengths,
            spectral_image,
            n_px, n_moments,
            32,
            n_bits,
            quantization_curve
        );
    }
}


// ----------------------------------------------------------------------------
// Curves alteration
// ----------------------------------------------------------------------------

void generate_random_curves(
    size_t lower_bound, size_t upper_bound,
    size_t n_curves_per_range,
    size_t n_elements,
    std::vector<std::vector<int>>& random_quantization_curves)
{
    random_quantization_curves.clear();
    random_quantization_curves.reserve((upper_bound - lower_bound + 1) * n_curves_per_range);

    for (size_t n_bits = lower_bound; n_bits <= upper_bound; n_bits++) {
        for (size_t iter = 0; iter < n_curves_per_range; iter++) {
            random_quantization_curves.resize(random_quantization_curves.size() + 1);
            std::vector<int>& quantization_curve = random_quantization_curves.back();
            quantization_curve.resize(n_elements);

            quantization_curve[0] = 32;

            for (size_t m = 1; m < n_elements; m++) {
                quantization_curve[m] = std::round((n_bits - 1) * (float)std::rand() / (float)RAND_MAX) + 1;
            }
        }
    }
}


void generate_mutated_curves(
    const std::vector<std::vector<int>>& original_curves,
    size_t n_mutations,
    std::vector<std::vector<int>>& mutated_curves)
{
    mutated_curves.clear();
    mutated_curves.reserve(original_curves.size() * n_mutations);

    for (const std::vector<int>& org_curve: original_curves) {
        for (size_t mutation = 0; mutation < n_mutations; mutation++) {
            mutated_curves.resize(mutated_curves.size() + 1);
            std::vector<int>& quantization_curve = mutated_curves.back();
            quantization_curve.resize(org_curve.size());

            quantization_curve[0] = org_curve[0];

            for (size_t m = 1; m < org_curve.size(); m++) {
                const int bits = std::round(16.f * ((float)std::rand() / (float)RAND_MAX) - 8.f);

                quantization_curve[m] = std::max(org_curve[m] + bits, 1);
            }
        }
    }
}


// ----------------------------------------------------------------------------
// Error computation for a set of quantization curves
// ----------------------------------------------------------------------------

void bounded_errors_for_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::vector<int>>& quantization_curves,
    std::vector<std::pair<int, double>>& errors)
{
    errors.clear();
    errors.reserve(quantization_curves.size());

    for (const std::vector<int>& quantization_curve: quantization_curves) {
        int n_bits = 32;

        for (size_t i = 1; i < quantization_curve.size(); i++) {
            n_bits += quantization_curve[i];
        }

        const double err = bounded_error_for_quantization_curve(
            wavelengths, spectral_image,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            quantization_curve);

        errors.push_back(std::make_pair(n_bits, err));
    }
}


void unbounded_errors_for_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::vector<int>>& quantization_curves,
    std::vector<std::pair<int, double>>& errors)
{
    errors.clear();
    errors.reserve(quantization_curves.size());

    for (const std::vector<int>& quantization_curve: quantization_curves) {
        int n_bits = 32;

        for (size_t i = 1; i < quantization_curve.size(); i++) {
            n_bits += quantization_curve[i];
        }

        const double err = unbounded_error_for_quantization_curve(
            wavelengths, spectral_image,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            quantization_curve);

        errors.push_back(std::make_pair(n_bits, err));
    }
}


void unbounded_to_bounded_errors_for_quantization_curves(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_px, size_t n_moments,
    const std::vector<double>& normalized_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    const std::vector<std::vector<int>>& quantization_curves,
    std::vector<std::pair<int, double>>& errors)
{
    errors.clear();
    errors.reserve(quantization_curves.size());

    for (const std::vector<int>& quantization_curve: quantization_curves) {
        int n_bits = 32;

        for (size_t i = 1; i < quantization_curve.size(); i++) {
            n_bits += quantization_curve[i];
        }

        const double err = unbounded_to_bounded_error_for_quantization_curve(
            wavelengths, spectral_image,
            n_px, n_moments,
            normalized_moments,
            mins, maxs,
            quantization_curve);

        errors.push_back(std::make_pair(n_bits, err));
    }
}



// ----------------------------------------------------------------------------
// Gnuplot output
// ----------------------------------------------------------------------------

void gnuplot_output_quantization_curves(
    const std::vector<std::vector<int>>& quantization_curves,
    std::ostream& output)
{
    for (const std::vector<int>& q_curve: quantization_curves) {
        for (size_t m = 1; m < q_curve.size(); m++) {
            output << m << " " << q_curve[m] << std::endl;
        }

        output << std::endl << std::endl;
    }
}


void gnuplot_output_quantization_curves(
    const std::vector<std::vector<int>>& quantization_curves,
    const char* filename)
{
    std::ofstream file(filename);

    if (file.is_open()) {
        gnuplot_output_quantization_curves(
            quantization_curves,
            file
        );

        file.close();
    } else {
        std::cerr << "Could not open " << filename << " for writing" << std::endl;
    }
}


void gnuplot_output_bits_errors(
    const std::vector<std::pair<int, double>> errors,
    std::ostream& output)
{
    for (const std::pair<int, double>& error: errors) {
        output << error.first << " " << error.second << std::endl;
    }
}


void gnuplot_output_bits_errors(
    const std::vector<std::pair<int, double>> errors_generated_curves,
    const std::vector<std::pair<int, double>> errors_mutated_curves,
    const std::vector<std::pair<int, double>> errors_random_curves,
    const char* filename)
{
    std::ofstream file(filename);

    if (file.is_open()) {
        gnuplot_output_bits_errors(errors_generated_curves, file);
        file << std::endl << std::endl;

        gnuplot_output_bits_errors(errors_mutated_curves, file);
        file << std::endl << std::endl;

        gnuplot_output_bits_errors(errors_random_curves, file);
        file << std::endl << std::endl;

        file.close();
    } else {
        std::cerr << "Could not open " << filename << " for writing" << std::endl;
    }
}


// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    (void)argc; (void)argv;

    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_spectral_image>" << std::endl;

        exit(0);
    }

    EXRSpectralImage image(argv[1]);

    const std::vector<SpectralFramebuffer*>& framebuffers = image.getSpectralFramebuffers();
    const SpectralFramebuffer* main_fb = framebuffers.front();

    std::cout << "Opened " << argv[1] << std::endl;
    std::cout << "Found: " << framebuffers.size() << " spectral framebuffer(s)." << std::endl;
    std::cout << "Using: " << main_fb->root_name << std::endl;

    const size_t n_wl = main_fb->wavelengths_nm.size();
    const size_t n_px = image.width() * image.height();
    const size_t n_moments = n_wl;

    std::vector<double> wavelengths(n_wl);
    for (size_t i = 0; i < n_wl; i++) {
        wavelengths[i] = main_fb->wavelengths_nm[i];
    }

    std::vector<double> spectral_image(n_wl * n_px);
    for (size_t wl = 0; wl < n_wl; wl++) {
        for (size_t px = 0; px < n_px; px++) {
            spectral_image[px * n_wl + wl] = main_fb->image_data[px * n_wl + wl];
        }
    }

    // Hard coded version

    // const size_t n_wl = 36;
    // const size_t n_px = 24;
    // const size_t n_moments = n_wl;

    // std::vector<double> wavelengths(n_wl);
    // for (size_t i = 0; i < n_wl; i++) {
    //     wavelengths[i] = macbeth_wavelengths[i];
    // }

    // std::vector<double> spectral_image(n_wl * n_px);
    // for (size_t wl = 0; wl < n_wl; wl++) {
    //     for (size_t px = 0; px < n_px; px++) {
    //         spectral_image[px * n_wl + wl] = macbeth_patches[px][wl];
    //     }
    // }

    // ------------------------------------------------------------------------
    // Error for quantizing a given moment
    // ------------------------------------------------------------------------

    // for (int m = 1; m < n_moments; m++) {
    //     quantize_dequantize_single_image(normalized_moments, n_px, n_moments, quantized_moments, m, 8);
    //     const double err = average_err(wavelengths, macbeth, n_px, n_moments, quantized_moments, mins, maxs);
    //     std::cout << m << " " << err << std::endl;
    // }

    // ------------------------------------------------------------------------
    // Determining quantization curves
    // ------------------------------------------------------------------------

    std::vector<std::vector<int>> bounded_quantization_curves;

    bounded_generate_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        bounded_quantization_curves
    );

    gnuplot_output_quantization_curves(
        bounded_quantization_curves,
        "q_curves_bounded.dat"
    );

    std::vector<std::vector<int>> unbounded_quantization_curves;

    unbounded_generate_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_quantization_curves
    );

    gnuplot_output_quantization_curves(
        unbounded_quantization_curves,
        "q_curves_unbounded.dat"
    );

    std::vector<std::vector<int>> unbounded_to_bounded_quantization_curves;

    unbounded_to_bounded_generate_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_to_bounded_quantization_curves
    );

    gnuplot_output_quantization_curves(
        unbounded_to_bounded_quantization_curves,
        "q_curves_unbounded_to_bounded.dat"
    );

    // ------------------------------------------------------------------------
    // Pareto
    // ------------------------------------------------------------------------

    std::vector<double> phases(n_wl);

    wavelengths_to_phases(wavelengths, phases);

    std::vector<double> moments(n_wl * n_px);

    compute_moments_image(
        phases, spectral_image,
        n_px,
        n_moments,
        moments
    );

    std::vector<double> bounded_compressed_moments(n_wl * n_px);
    std::vector<double> bounded_normalized_moments(n_wl * n_px);
    std::vector<double> bounded_mins, bounded_maxs;

    bounded_compress_moments_image(
        moments,
        n_px,
        n_moments,
        bounded_compressed_moments
    );

    normalize_moment_image(
        bounded_compressed_moments,
        n_px, n_moments,
        bounded_normalized_moments,
        bounded_mins, bounded_maxs
    );

    std::vector<double> unbounded_compressed_moments(n_wl * n_px);
    std::vector<double> unbounded_normalized_moments(n_wl * n_px);
    std::vector<double> unbounded_mins, unbounded_maxs;

    unbounded_compress_moments_image(
        moments,
        n_px,
        n_moments,
        unbounded_compressed_moments
    );

    normalize_moment_image(
        unbounded_compressed_moments,
        n_px,
        n_moments,
        unbounded_normalized_moments,
        unbounded_mins, unbounded_maxs
    );

    std::vector<double> unbounded_to_bounded_compressed_moments(n_wl * n_px);
    std::vector<double> unbounded_to_bounded_normalized_moments(n_wl * n_px);
    std::vector<double> unbounded_to_bounded_mins, unbounded_to_bounded_maxs;

    unbounded_to_bounded_compress_moments_image(
        moments,
        n_px,
        n_moments,
        unbounded_to_bounded_compressed_moments
    );

    normalize_moment_image(
        unbounded_to_bounded_compressed_moments,
        n_px,
        n_moments,
        unbounded_to_bounded_normalized_moments,
        unbounded_to_bounded_mins, unbounded_to_bounded_maxs
    );

    // Bounded
    std::vector<std::vector<int>> bounded_mutated_curves;
    std::vector<std::vector<int>> bounded_random_curves;

    generate_mutated_curves(bounded_quantization_curves, 500, bounded_mutated_curves);
    generate_random_curves(3, 16, 500, n_moments, bounded_random_curves);

    std::vector<std::pair<int, double>> error_bounded_quantization_curves;
    std::vector<std::pair<int, double>> error_bounded_mutated_curves;
    std::vector<std::pair<int, double>> error_bounded_random_curves;

    bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        bounded_normalized_moments, bounded_mins, bounded_maxs,
        bounded_quantization_curves,
        error_bounded_quantization_curves
    );

    bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        bounded_normalized_moments, bounded_mins, bounded_maxs,
        bounded_mutated_curves,
        error_bounded_mutated_curves
    );

    bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        bounded_normalized_moments, bounded_mins, bounded_maxs,
        bounded_random_curves,
        error_bounded_random_curves
    );

    gnuplot_output_bits_errors(
        error_bounded_quantization_curves,
        error_bounded_mutated_curves,
        error_bounded_random_curves,
        "pareto_bounded.dat"
    );

    // Unbounded
    std::vector<std::vector<int>> unbounded_mutated_curves;
    std::vector<std::vector<int>> unbounded_random_curves;

    generate_mutated_curves(unbounded_quantization_curves, 500, unbounded_mutated_curves);
    generate_random_curves(3, 16, 500, n_moments, unbounded_random_curves);

    std::vector<std::pair<int, double>> error_unbounded_quantization_curves;
    std::vector<std::pair<int, double>> error_unbounded_mutated_curves;
    std::vector<std::pair<int, double>> error_unbounded_random_curves;

    unbounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_normalized_moments, unbounded_mins, unbounded_maxs,
        unbounded_quantization_curves,
        error_unbounded_quantization_curves
    );

    unbounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_normalized_moments, unbounded_mins, unbounded_maxs,
        unbounded_mutated_curves,
        error_unbounded_mutated_curves
    );

    unbounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_normalized_moments, unbounded_mins, unbounded_maxs,
        unbounded_random_curves,
        error_unbounded_random_curves
    );

    gnuplot_output_bits_errors(
        error_unbounded_quantization_curves,
        error_unbounded_mutated_curves,
        error_unbounded_random_curves,
        "pareto_unbounded.dat"
    );

    // Unbounded to bounded
    std::vector<std::vector<int>> unbounded_to_bounded_mutated_curves;
    std::vector<std::vector<int>> unbounded_to_bounded_random_curves;

    generate_mutated_curves(unbounded_to_bounded_quantization_curves, 500, unbounded_to_bounded_mutated_curves);
    generate_random_curves(3, 16, 500, n_moments, unbounded_to_bounded_random_curves);

    std::vector<std::pair<int, double>> error_unbounded_to_bounded_quantization_curves;
    std::vector<std::pair<int, double>> error_unbounded_to_bounded_mutated_curves;
    std::vector<std::pair<int, double>> error_unbounded_to_bounded_random_curves;

    unbounded_to_bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_to_bounded_normalized_moments, unbounded_to_bounded_mins, unbounded_to_bounded_maxs,
        unbounded_to_bounded_quantization_curves,
        error_unbounded_to_bounded_quantization_curves
    );

    unbounded_to_bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_to_bounded_normalized_moments, unbounded_to_bounded_mins, unbounded_to_bounded_maxs,
        unbounded_to_bounded_mutated_curves,
        error_unbounded_to_bounded_mutated_curves
    );

    unbounded_to_bounded_errors_for_quantization_curves(
        wavelengths, spectral_image,
        n_px, n_moments,
        unbounded_to_bounded_normalized_moments, unbounded_to_bounded_mins, unbounded_to_bounded_maxs,
        unbounded_to_bounded_random_curves,
        error_unbounded_to_bounded_random_curves
    );

    gnuplot_output_bits_errors(
        error_unbounded_to_bounded_quantization_curves,
        error_unbounded_to_bounded_mutated_curves,
        error_unbounded_to_bounded_random_curves,
        "pareto_unbounded_to_bounded.dat"
    );

    return 0;
}

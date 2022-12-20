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

#include <EXRSpectralImage.h>
#include <curve_quantization.h>
#include <sstream>
#include <fstream>
#include <chrono>

#define VERBOSE

void run_for_bounded(
    const std::vector<double>& wavelengths,
    const std::vector<double>& image_data,
    int n_pixels, int n_bands,
    int n_bits_dc,
    int n_bits_ac1,
    const std::string& output_prefix)
{
#ifdef VERBOSE
    std::cout << "Running for bounded...";
    auto start = std::chrono::steady_clock::now();
#endif // VERBOSE

    std::vector<int> quantization_curve_b;

    double err_utb = bounded_compute_quantization_curve(
        wavelengths,
        image_data,
        n_pixels, n_bands,
        n_bits_dc,
        n_bits_ac1,
        quantization_curve_b
    );

    // Save data
    std::stringstream output_file;
    output_file << output_prefix << "_b_" << n_bits_ac1 << ".txt";

    std::ofstream out_b(output_file.str());

    for (size_t i = 0; i < quantization_curve_b.size(); i++) {
        out_b << quantization_curve_b[i] << " ";
    }

    out_b << std::endl << err_utb;

#ifdef VERBOSE
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "\t\t\t" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    std::cout << "Writing: " << output_file.str() << std::endl;
#endif // VERBOSE
}


void run_for_unbounded(
    const std::vector<double>& wavelengths,
    const std::vector<double>& image_data,
    int n_pixels, int n_bands,
    int n_bits_dc,
    int n_bits_ac1,
    const std::string& output_prefix)
{
#ifdef VERBOSE
    std::cout << "Running for unbounded...";
    auto start = std::chrono::steady_clock::now();
#endif // VERBOSE

    std::vector<int> quantization_curve_u;

    double err_utb = unbounded_compute_quantization_curve(
        wavelengths,
        image_data,
        n_pixels, n_bands,
        n_bits_dc,
        n_bits_ac1,
        quantization_curve_u
    );

    // Save data
    std::stringstream output_file;
    output_file << output_prefix << "_u_" << n_bits_ac1 << ".txt";

    std::ofstream out_u(output_file.str());

    for (size_t i = 0; i < quantization_curve_u.size(); i++) {
        out_u << quantization_curve_u[i] << " ";
    }

    out_u << std::endl << err_utb;

#ifdef VERBOSE
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "\t\t" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    std::cout << "Writing: " << output_file.str() << std::endl;
#endif // VERBOSE
}


void run_for_unbounded_to_bounded(
    const std::vector<double>& wavelengths,
    const std::vector<double>& image_data,
    int n_pixels, int n_bands,
    int n_bits_dc,
    int n_bits_ac1,
    const std::string& output_prefix)
{
#ifdef VERBOSE
    std::cout << "Running for unbounded to bounded...";
    auto start = std::chrono::steady_clock::now();
#endif // VERBOSE

    std::vector<int> quantization_curve_utb;

    double err_utb = unbounded_to_bounded_compute_quantization_curve(
        wavelengths,
        image_data,
        n_pixels, n_bands,
        n_bits_dc,
        n_bits_ac1,
        quantization_curve_utb
    );

    // Save data
    std::stringstream output_file;
    output_file << output_prefix << "_utb_" << n_bits_ac1 << ".txt";

    std::ofstream out_utb(output_file.str());

    for (size_t i = 0; i < quantization_curve_utb.size(); i++) {
        out_utb << quantization_curve_utb[i] << " ";
    }

    out_utb << std::endl << err_utb;

#ifdef VERBOSE
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "\t" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    std::cout << "Writing: " << output_file.str() << std::endl;
#endif // VERBOSE
}


void run_for_upperbound(
    const std::vector<double>& wavelengths,
    const std::vector<double>& image_data,
    int n_pixels, int n_bands,
    int n_bits_dc,
    int n_bits_ac1,
    const std::string& output_prefix)
{
#ifdef VERBOSE
    std::cout << "Running for upperbound...";
    auto start = std::chrono::steady_clock::now();
#endif // VERBOSE

    std::vector<int> quantization_curve_utb;

    double err_utb = upperbound_compute_quantization_curve(
        wavelengths,
        image_data,
        n_pixels, n_bands,
        n_bits_dc,
        n_bits_ac1,
        quantization_curve_utb
    );

    // Save data
    std::stringstream output_file;
    output_file << output_prefix << "_up_" << n_bits_ac1 << ".txt";

    std::ofstream out_utb(output_file.str());

    for (size_t i = 0; i < quantization_curve_utb.size(); i++) {
        out_utb << quantization_curve_utb[i] << " ";
    }

    // Add the relative scale memory footprint
    out_utb << "8 ";

    out_utb << std::endl << err_utb;

#ifdef VERBOSE
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "\t\t" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    std::cout << "Writing: " << output_file.str() << std::endl;
#endif // VERBOSE
}


void run_for_twobounds(
    const std::vector<double>& wavelengths,
    const std::vector<double>& image_data,
    int n_pixels, int n_bands,
    int n_bits_dc,
    int n_bits_ac1,
    const std::string& output_prefix)
{
#ifdef VERBOSE
    std::cout << "Running for twobounds...";
    auto start = std::chrono::steady_clock::now();
#endif // VERBOSE

    std::vector<int> quantization_curve_utb;

    double err_utb = twobounds_compute_quantization_curve(
        wavelengths,
        image_data,
        n_pixels, n_bands,
        n_bits_dc,
        n_bits_ac1,
        quantization_curve_utb
    );

    // Save data
    std::stringstream output_file;
    output_file << output_prefix << "_tb_" << n_bits_ac1 << ".txt";

    std::ofstream out_utb(output_file.str());

    for (size_t i = 0; i < quantization_curve_utb.size(); i++) {
        out_utb << quantization_curve_utb[i] << " ";
    }

    // Add the relative scale memory footprint
    out_utb << "8 ";

    out_utb << std::endl << err_utb;

#ifdef VERBOSE
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "\t\t" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    std::cout << "Writing: " << output_file.str() << std::endl;
#endif // VERBOSE
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral exr> <output_prefix> <run_for_bounded y/n>" << std::endl;
        return 0;
    }

    const char* filename = argv[1];
    const std::string output_prefix = argv[2];

#ifdef VERBOSE
    std::cout << "Benchmarking: " << filename << std::endl;
#endif // VERBOSE

    // Optional run: not relevant for spectra with a value > 1
    bool do_run_for_bounded = true;

    if (argc == 4) {
        if (argv[3][0] == 'n') {
            do_run_for_bounded = false;
        }
    }

    EXRSpectralImage image(filename);

    const int n_pixels = image.width() * image.height();

    // const int n_bits[] = {12, 10, 8};
    const int n_bits[] = {6, 7, 8, 9, 10, 11, 12, 13};

    for (SpectralFramebuffer* fb: image.getSpectralFramebuffers()) {
        for (int n_bits_ac1: n_bits) {
            #ifdef VERBOSE
            std::cout << "Using " << n_bits_ac1 << "bits" << std::endl;
            #endif // VERBOSE

            const int n_bands = fb->wavelengths_nm.size();

            std::vector<double> wavelengths(n_bands);
            std::vector<double> image_data(n_pixels * n_bands);

            // Cast to double
            for (int i = 0; i < n_bands; i++) {
                wavelengths[i] = (double)fb->wavelengths_nm[i];
            }

            for (int i = 0; i < n_pixels * n_bands; i++) {
                image_data[i] = (double)fb->image_data[i];
            }

            int n_bits_dc = fb->pixel_type == PixelType::HALF ? 16 : 32;

            if (do_run_for_bounded) {
                run_for_bounded(
                    wavelengths,
                    image_data,
                    n_pixels, n_bands,
                    n_bits_dc,
                    n_bits_ac1,
                    output_prefix
                );
            }

            run_for_unbounded(
                wavelengths,
                image_data,
                n_pixels, n_bands,
                n_bits_dc,
                n_bits_ac1,
                output_prefix
            );

            run_for_unbounded_to_bounded(
                wavelengths,
                image_data,
                n_pixels, n_bands,
                n_bits_dc,
                n_bits_ac1,
                output_prefix
            );

            run_for_upperbound(
                wavelengths,
                image_data,
                n_pixels, n_bands,
                n_bits_dc,
                n_bits_ac1,
                output_prefix
            );

            run_for_twobounds(
                wavelengths,
                image_data,
                n_pixels, n_bands,
                n_bits_dc,
                n_bits_ac1,
                output_prefix
            );
        }
    }

    return 0;
}

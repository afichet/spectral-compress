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
#include <JXLImage.h>

#include <moments_image.h>


void decompress_spectral_framebuffer(
    const std::vector<float>& wavelengths,
    const std::vector<std::vector<float>>& compressed_moments,
    const std::vector<float>& mins,
    const std::vector<float>& maxs,
    const std::vector<uint8_t>& relative_scales,
    float global_max,
    std::vector<float>& spectral_framebuffer)
{
    const size_t n_moments = compressed_moments.size();
    const size_t n_pixels = compressed_moments[0].size();

    std::vector<double> wavelengths_d(wavelengths.size());
    std::vector<double> compressed_moments_d(n_moments * n_pixels);
    std::vector<double> mins_d(n_moments - 1), maxs_d(n_moments - 1);
    std::vector<double> spectral_framebuffer_d;

    for (size_t i = 0; i < wavelengths.size(); i++) {
        wavelengths_d[i] = wavelengths[i];
    }

    for (size_t px = 0; px < n_pixels; px++) {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments_d[px * n_moments + m] = compressed_moments[m][px];
        }
    }

    for (size_t i = 0; i < n_moments - 1; i++) {
        mins_d[i] = mins[i];
        maxs_d[i] = maxs[i];
    }

    upperbound_decompress_spectral_image(
        wavelengths_d,
        compressed_moments_d,
        mins_d, maxs_d,
        relative_scales,
        global_max,
        n_pixels, n_moments,
        spectral_framebuffer_d
    );

    spectral_framebuffer.resize(spectral_framebuffer_d.size());

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_framebuffer_d.size(); i++) {
        spectral_framebuffer[i] = spectral_framebuffer_d[i];
    }
}

template<typename T>
T error_images(
    const std::vector<T>& reference,
    const std::vector<T>& comparison,
    size_t n_pixels,
    size_t n_bands)
{
    T error = 0;

    for (size_t p = 0; p < n_pixels; p++) {
        T px_sum_err = 0;
        T avg = 0;

        for (size_t i = 0; i < n_bands; i++) {
            const T q = reference[p * n_bands + i] - comparison[p * n_bands + i];
            avg += reference[p * n_bands + i];
            px_sum_err += q * q;
        }

        avg /= n_bands;

        if (avg > 0) {
            error += std::sqrt(px_sum_err) / avg;
        }
    }

    return error / (double)n_pixels;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << "<ref_dump> <compressed_dump>" << std::endl;
        return 0;
    }

    // Load both dumps
    const char* ref_dump = argv[1];
    const char* cmp_dump = argv[2];


    // std::string ref_dump_str = argv[1];
    // std::string cmp_dump_str = argv[2];

    // if (ref_dump_str.)

    EXRSpectralImage* exr_image = EXRSpectralImage::read_dump(ref_dump);
    JXLImage*         jxl_image = JXLImage::read_dump(cmp_dump);

    // Retrieve original spectral framebuffer
    SpectralFramebuffer* fb = exr_image->getSpectralFramebuffers()[0];
    const std::vector<float>& spectral_framebuffer_original = fb->image_data;

    // We need to decompress cmp
    const SGEGBox box   = jxl_image->getBox();
    const size_t width  = jxl_image->width();
    const size_t height = jxl_image->height();
    const size_t n_bands = fb->wavelengths_nm.size();

    SGEGSpectralGroup sg = box.spectral_groups[0];
    std::string root_name = sg.root_name.data();
    // 1 extra framebuffer is used to store relative max
    const size_t n_moments = sg.layer_indices.size() - 1;

    std::vector<std::vector<float>> moments(n_moments);
    std::vector<float> relative_scales_f;
    std::vector<float> spectral_framebuffer_decompressed;

    for (size_t m = 0; m < n_moments; m++) {
        moments[m] = jxl_image->getFramebufferDataConst(sg.layer_indices[m]);
    }

    relative_scales_f = jxl_image->getFramebufferDataConst(sg.layer_indices.back());

    // TODO: This is a bit hacky
    std::vector<uint8_t> relative_scales(relative_scales_f.size());

    for (size_t i = 0; i < relative_scales_f.size(); i++) {
        relative_scales[i] = std::numeric_limits<uint8_t>::max() * relative_scales_f[i];
    }

    decompress_spectral_framebuffer(
        sg.wavelengths,
        moments,
        sg.mins, sg.maxs,
        relative_scales,
        sg.global_max,
        spectral_framebuffer_decompressed
    );


    // Compute the error
    double err = 0;

    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            double curr_err = 0;
            double avg = 0;

            for (size_t band = 0; band < n_bands; band++) {
                float ref = spectral_framebuffer_original[(y * width + x) * n_bands + band];
                float cmp = spectral_framebuffer_decompressed[(y * width + x) * n_bands + band];

                curr_err += (ref - cmp) * (ref - cmp);

                avg += ref;
            }

            if (avg > 0) {
                err += std::sqrt(curr_err) / (avg / (float)n_bands);
            }
        }
    }


    err = err / (float)(width * height);

    std::cout << "Error: " << err << std::endl;


    const float err2 = error_images(
        spectral_framebuffer_original,
        spectral_framebuffer_decompressed,
        width * height, n_bands
    );

    std::cout << "Error2: " << err2 << std::endl;

    delete exr_image;
    delete jxl_image;

    return 0;
}

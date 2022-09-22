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

#include "moments_image.h"
#include "moments.h"

#include <cstring>
#include <Eigen/Core>

/******************************************************************************
 * C interface
 *****************************************************************************/

extern "C" {

void compute_moments_image(
    const float phases[],
    size_t n_phases,
    const float spectral_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[])
{
    // TODO: Simplify
    if (n_moments == n_phases) {
        std::vector<float> transform(n_phases * n_phases);
        compute_basis_signal_to_moments(phases, n_phases, transform.data());

        Eigen::Map<Eigen::MatrixXf> t_mat(transform.data(), n_phases, n_phases);

        // Local copy for non const manipulation
        std::vector<float> spec_image(width * height * n_phases);
        std::memcpy(spec_image.data(), spectral_image, width * height * n_phases * sizeof(float));

        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            const Eigen::Map<Eigen::VectorXf> signal(&(spec_image[n_phases * i]), n_phases);
            Eigen::Map<Eigen::VectorXf> moments(&(moments_image[n_phases * i]), n_phases);

            moments = t_mat * signal;
        }
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            compute_moments(
                phases, 
                n_phases, 
                &(spectral_image[n_phases * i]),
                n_moments,
                &(moments_image[(n_moments + 1) * i]));
        }
    }
}


void compress_moments_image(
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float compressed_moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        compress_moments(
            &(moments_image[n_moments * i]),
            n_moments,
            &(compressed_moments_image[n_moments * i])
        );
    }
}


void compress_bounded_moments_image(
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float compressed_moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        compress_bounded_moments(
            &(moments_image[n_moments * i]),
            n_moments,
            &(compressed_moments_image[n_moments * i])
        );
    }
}


void decompress_moments_image(
    const float compressed_moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        decompress_moments(
            &(compressed_moments_image[n_moments * i]),
            n_moments,
            &(moments_image[n_moments * i])
        );
    }
}


void decompress_bounded_moments_image(
    const float compressed_moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        decompress_bounded_moments(
            &(compressed_moments_image[n_moments * i]),
            n_moments,
            &(moments_image[n_moments * i])
        );
    }
}


void compute_density_image(
    const float phases[],
    size_t n_phases,
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float density_image[])
{
    // TODO: Simplify
    if (n_moments == n_phases) {
        std::vector<float> transform(n_phases * n_phases);
        compute_basis_moments_to_signal(phases, n_phases, transform.data());

        Eigen::Map<Eigen::MatrixXf> t_mat(transform.data(), n_phases, n_phases);

        // Local copy for non const manipulation
        std::vector<float> m_img(width * height * n_moments);
        std::memcpy(m_img.data(), moments_image, width * height * n_moments * sizeof(float));

        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            const Eigen::Map<Eigen::VectorXf> moments(&(m_img[n_phases * i]), n_phases);
            Eigen::Map<Eigen::VectorXf> signal(&(density_image[n_phases * i]), n_phases);

            signal = t_mat * moments;
        }
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < width * height; i++) {
            compute_density(
                phases,
                n_phases,
                &(moments_image[n_moments * i]),
                n_moments,
                &(density_image[n_phases * i])
            );
        }
    }
}


void compute_density_bounded_lagrange_image(
    const float phases[],
    size_t n_phases,
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float density_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        compute_density_bounded_lagrange(
            phases,
            n_phases,
            &(moments_image[n_moments * i]),
            n_moments,
            &(density_image[n_phases * i])
        );
    }
}


void normalize_moment_image(
    const float src[],
    size_t n_px, size_t n_moments,
    float dest[],
    float mins[],
    float maxs[])
{    
    // Initialize vectors
    for (size_t m = 1; m < n_moments; m++) {
        mins[m - 1] = src[m];
        maxs[m - 1] = src[m];
    }

    // Find mins & maxs
    for (size_t px = 0; px < n_px; px++) {
        for (size_t m = 1; m < n_moments; m++) {
            mins[m - 1] = std::min(mins[m - 1], src[px * n_moments + m]);
            maxs[m - 1] = std::max(maxs[m - 1], src[px * n_moments + m]);
        }
    }

    // Rescale
    #pragma omp parallel for
    for (size_t px = 0; px < n_px; px++) {
        // DC component
        dest[px * n_moments] = src[px * n_moments];

        // AC components
        for (size_t m = 1; m < n_moments; m++) {
            const float& min = mins[m - 1];
            const float& max = maxs[m - 1];

            dest[px * n_moments + m] = (src[px * n_moments + m] - min) / (max - min);
        }
    }
}


void denormalize_moment_image(
    const float src[],
    size_t n_px, size_t n_moments,
    const float mins[],
    const float maxs[],
    float dest[])
{
    #pragma omp parallel for
    for (size_t px = 0; px < n_px; px++) {
        // DC component
        dest[px * n_moments] = src[px * n_moments];

        // AC components
        for (size_t m = 1; m < n_moments; m++) {
            const float& min = mins[m - 1];
            const float& max = maxs[m - 1];

            dest[px * n_moments + m] = (max - min) * src[px * n_moments + m] + min;
        }
    }
}

} // extern "C"


/******************************************************************************
 * C++ interface
 *****************************************************************************/

void compute_moments_image(
    const std::vector<float>& phases,
    const std::vector<float>& spectral_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& moments_image)
{
    moments_image.resize(width * height * n_moments);

    compute_moments_image(
        phases.data(),
        phases.size(),
        spectral_image.data(),
        width, height,
        n_moments,
        moments_image.data());
}


void compress_moments_image(
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& compressed_moments_image)
{
    compressed_moments_image.resize(moments_image.size());

    compress_moments_image(
        moments_image.data(), 
        width, height, 
        n_moments, 
        compressed_moments_image.data()
    );
}


void compress_bounded_moments_image(
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& compressed_moments_image)
{
    compressed_moments_image.resize(moments_image.size());

    compress_bounded_moments_image(
        moments_image.data(), 
        width, height, 
        n_moments, 
        compressed_moments_image.data()
    );
}


void decompress_moments_image(
    const std::vector<float>& compressed_moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& moments_image)
{
    moments_image.resize(compressed_moments_image.size());

    decompress_moments_image(
        compressed_moments_image.data(),
        width, height,
        n_moments,
        moments_image.data()
    );
}


void decompress_bounded_moments_image(
    const std::vector<float>& compressed_moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& moments_image)
{
    moments_image.resize(compressed_moments_image.size());

    decompress_bounded_moments_image(
        compressed_moments_image.data(),
        width, height,
        n_moments,
        moments_image.data()
    );
}


void compute_density_image(
    const std::vector<float>& phases,
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& density_image)
{
    density_image.resize(phases.size() * width * height);

    compute_density_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        width, height,
        n_moments,
        density_image.data()
    );
}


void compute_density_bounded_lagrange_image(
    const std::vector<float>& phases,
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& density_image)
{
    density_image.resize(phases.size() * width * height);

    compute_density_bounded_lagrange_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        width, height,
        n_moments,
        density_image.data()
    );
}


void normalize_moment_image(
    const std::vector<float>& src,
    size_t n_px, size_t n_moments,
    std::vector<float>& dest,
    std::vector<float>& mins,
    std::vector<float>& maxs)
{
    dest.resize(src.size());
    mins.resize(n_px * (n_moments - 1));
    maxs.resize(n_px * (n_moments - 1));

    normalize_moment_image(
        src.data(),
        n_px, n_moments,
        dest.data(),
        mins.data(),
        maxs.data()
    );
}


void denormalize_moment_image(
    const std::vector<float>& src,
    size_t n_px, size_t n_moments,
    const std::vector<float>& mins,
    const std::vector<float>& maxs,
    std::vector<float>& dest)
{
    dest.resize(src.size());

    denormalize_moment_image(
        src.data(),
        n_px, n_moments,
        mins.data(),
        maxs.data(),
        dest.data()
    );
}
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


extern "C"
void compute_moments_image(
    const double phases[],
    size_t n_phases,
    const double spectral_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[])
{
    // TODO: Simplify
    if (n_moments == n_phases) {
        std::vector<double> transform(n_phases * n_phases);
        compute_basis_signal_to_moments(phases, n_phases, transform.data());

        Eigen::Map<Eigen::MatrixXd> t_mat(transform.data(), n_phases, n_phases);

        // Local copy for non const manipulation
        std::vector<double> spec_image(n_pixels * n_phases);
        std::memcpy(spec_image.data(), spectral_image, n_pixels * n_phases * sizeof(double));

        #pragma omp parallel for
        for (size_t i = 0; i < n_pixels; i++) {
            const Eigen::Map<Eigen::VectorXd> signal(&(spec_image[n_phases * i]), n_phases);
            Eigen::Map<Eigen::VectorXd> moments(&(moments_image[n_phases * i]), n_phases);

            moments = t_mat * signal;
        }
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < n_pixels; i++) {
            compute_moments(
                phases,
                n_phases,
                &(spectral_image[n_phases * i]),
                n_moments,
                &(moments_image[(n_moments + 1) * i]));
        }
    }
}

void compute_moments_image(
    const std::vector<double>& phases,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image)
{
    moments_image.resize(n_pixels * n_moments);

    compute_moments_image(
        phases.data(),
        phases.size(),
        spectral_image.data(),
        n_pixels,
        n_moments,
        moments_image.data());
}


/*****************************************************************************/
/* Compression                                                               */
/*****************************************************************************/

extern "C"
void unbounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        unbounded_compress_moments(
            &(moments_image[n_moments * i]),
            n_moments,
            &(compressed_moments_image[n_moments * i])
        );
    }
}

void unbounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compressed_moments_image)
{
    compressed_moments_image.resize(moments_image.size());

    unbounded_compress_moments_image(
        moments_image.data(),
        n_pixels,
        n_moments,
        compressed_moments_image.data()
    );
}


extern "C"
void bounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        bounded_compress_moments(
            &(moments_image[n_moments * i]),
            n_moments,
            &(compressed_moments_image[n_moments * i])
        );
    }
}

void bounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compressed_moments_image)
{
    compressed_moments_image.resize(moments_image.size());

    bounded_compress_moments_image(
        moments_image.data(),
        n_pixels,
        n_moments,
        compressed_moments_image.data()
    );
}


extern "C"
void unbounded_to_bounded_compress_moments_image(
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double compressed_moment_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        unbounded_to_bounded_compress_moments(
            &(moments_image[n_moments * i]),
            n_moments,
            &(compressed_moment_image[n_moments * i])
        );
    }
}

void unbounded_to_bounded_compress_moments_image(
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& compressed_moments_image)
{
    compressed_moments_image.resize(moments_image.size());

    unbounded_to_bounded_compress_moments_image(
        moments_image.data(),
        n_pixels,
        n_moments,
        compressed_moments_image.data()
    );
}


/*****************************************************************************/
/* Decompression                                                             */
/*****************************************************************************/

extern "C"
void unbounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        unbounded_decompress_moments(
            &(compressed_moments_image[n_moments * i]),
            n_moments,
            &(moments_image[n_moments * i])
        );
    }
}

void unbounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image)
{
    moments_image.resize(compressed_moments_image.size());

    unbounded_decompress_moments_image(
        compressed_moments_image.data(),
        n_pixels,
        n_moments,
        moments_image.data()
    );
}


extern "C"
void bounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        bounded_decompress_moments(
            &(compressed_moments_image[n_moments * i]),
            n_moments,
            &(moments_image[n_moments * i])
        );
    }
}

void bounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image)
{
    moments_image.resize(compressed_moments_image.size());

    bounded_decompress_moments_image(
        compressed_moments_image.data(),
        n_pixels,
        n_moments,
        moments_image.data()
    );
}


extern "C"
void unbounded_to_bounded_decompress_moments_image(
    const double compressed_moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double moments_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        unbounded_to_bounded_decompress_moments(
            &(compressed_moments_image[n_moments * i]),
            n_moments,
            &(moments_image[n_moments * i])
        );
    }
}

void unbounded_to_bounded_decompress_moments_image(
    const std::vector<double>& compressed_moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& moments_image)
{
    moments_image.resize(compressed_moments_image.size());

    unbounded_to_bounded_decompress_moments_image(
        compressed_moments_image.data(),
        n_pixels,
        n_moments,
        moments_image.data()
    );
}


/*****************************************************************************/

extern "C"
void compute_density_image(
    const double phases[],
    size_t n_phases,
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double density_image[])
{
    // TODO: Simplify
    if (n_moments == n_phases) {
        std::vector<double> transform(n_phases * n_phases);
        compute_basis_moments_to_signal(phases, n_phases, transform.data());

        Eigen::Map<Eigen::MatrixXd> t_mat(transform.data(), n_phases, n_phases);

        // Local copy for non const manipulation
        std::vector<double> m_img(n_pixels * n_moments);
        std::memcpy(m_img.data(), moments_image, n_pixels * n_moments * sizeof(double));

        #pragma omp parallel for
        for (size_t i = 0; i < n_pixels; i++) {
            const Eigen::Map<Eigen::VectorXd> moments(&(m_img[n_phases * i]), n_phases);
            Eigen::Map<Eigen::VectorXd> signal(&(density_image[n_phases * i]), n_phases);

            signal = t_mat * moments;
        }
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < n_pixels; i++) {
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

void compute_density_image(
    const std::vector<double>& phases,
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& density_image)
{
    density_image.resize(phases.size() * n_pixels);

    compute_density_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        n_pixels,
        n_moments,
        density_image.data()
    );
}


extern "C"
void bounded_compute_density_lagrange_image(
    const double phases[],
    size_t n_phases,
    const double moments_image[],
    size_t n_pixels,
    size_t n_moments,
    double density_image[])
{
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        bounded_compute_density_lagrange(
            phases,
            n_phases,
            &(moments_image[n_moments * i]),
            n_moments,
            &(density_image[n_phases * i])
        );
    }
}

void bounded_compute_density_lagrange_image(
    const std::vector<double>& phases,
    const std::vector<double>& moments_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& density_image)
{
    density_image.resize(phases.size() * n_pixels);

    bounded_compute_density_lagrange_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        n_pixels,
        n_moments,
        density_image.data()
    );
}


extern "C"
void normalize_moment_image(
    const double src[],
    size_t n_pixels, size_t n_moments,
    double dest[],
    double mins[],
    double maxs[])
{
    // Initialize vectors
    for (size_t m = 1; m < n_moments; m++) {
        mins[m - 1] = src[m];
        maxs[m - 1] = src[m];
    }

    // Find mins & maxs
    for (size_t px = 0; px < n_pixels; px++) {
        for (size_t m = 1; m < n_moments; m++) {
            mins[m - 1] = std::min(mins[m - 1], src[px * n_moments + m]);
            maxs[m - 1] = std::max(maxs[m - 1], src[px * n_moments + m]);
        }
    }

    // Rescale
    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // DC component
        dest[px * n_moments] = src[px * n_moments];

        // AC components
        for (size_t m = 1; m < n_moments; m++) {
            const double& min = mins[m - 1];
            const double& max = maxs[m - 1];

            dest[px * n_moments + m] = (src[px * n_moments + m] - min) / (max - min);
        }
    }
}

void normalize_moment_image(
    const std::vector<double>& src,
    size_t n_px, size_t n_moments,
    std::vector<double>& dest,
    std::vector<double>& mins,
    std::vector<double>& maxs)
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


extern "C"
void denormalize_moment_image(
    const double src[],
    size_t n_pixels, size_t n_moments,
    const double mins[],
    const double maxs[],
    double dest[])
{
    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // DC component
        dest[px * n_moments] = src[px * n_moments];

        // AC components
        for (size_t m = 1; m < n_moments; m++) {
            const double& min = mins[m - 1];
            const double& max = maxs[m - 1];

            dest[px * n_moments + m] = (max - min) * src[px * n_moments + m] + min;
        }
    }
}

void denormalize_moment_image(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    std::vector<double>& dest)
{
    dest.resize(src.size());

    denormalize_moment_image(
        src.data(),
        n_pixels, n_moments,
        mins.data(),
        maxs.data(),
        dest.data()
    );
}




/*****************************************************************************/
/* Full pipeline                                                             */
/*****************************************************************************/

// Compression

void linear_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs)
{
    std::vector<double> phases;
    std::vector<double> moments_image;
    std::vector<double> compressed_moments_image;

    wavelengths_to_phases(wavelengths, phases);

    compute_moments_image(
        phases,
        spectral_image,
        n_pixels, n_moments,
        moments_image
    );

    normalize_moment_image(
        moments_image,
        n_pixels, n_moments,
        normalized_moments_image,
        mins, maxs
    );
}


void bounded_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs)
{
    std::vector<double> phases;
    std::vector<double> moments_image;
    std::vector<double> compressed_moments_image;

    wavelengths_to_phases(wavelengths, phases);

    compute_moments_image(
        phases,
        spectral_image,
        n_pixels, n_moments,
        moments_image
    );

    bounded_compress_moments_image(
        moments_image,
        n_pixels, n_moments,
        compressed_moments_image
    );

    normalize_moment_image(
        compressed_moments_image,
        n_pixels, n_moments,
        normalized_moments_image,
        mins, maxs
    );
}


void unbounded_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs)
{
    std::vector<double> phases;
    std::vector<double> moments_image;
    std::vector<double> compressed_moments_image;

    wavelengths_to_phases(wavelengths, phases);

    compute_moments_image(
        phases,
        spectral_image,
        n_pixels, n_moments,
        moments_image
    );

    unbounded_compress_moments_image(
        moments_image,
        n_pixels, n_moments,
        compressed_moments_image
    );

    normalize_moment_image(
        compressed_moments_image,
        n_pixels, n_moments,
        normalized_moments_image,
        mins, maxs
    );
}


void unbounded_to_bounded_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs)
{
    std::vector<double> phases;
    std::vector<double> moments_image;
    std::vector<double> compressed_moments_image;

    wavelengths_to_phases(wavelengths, phases);

    compute_moments_image(
        phases,
        spectral_image,
        n_pixels, n_moments,
        moments_image
    );

    unbounded_to_bounded_compress_moments_image(
        moments_image,
        n_pixels, n_moments,
        compressed_moments_image
    );

    normalize_moment_image(
        compressed_moments_image,
        n_pixels, n_moments,
        normalized_moments_image,
        mins, maxs
    );
}


// Decompression

void linear_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image)
{
    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image(
        normalized_moments_image,
        n_pixels,
        n_moments,
        mins, maxs,
        moments_image
    );

    compute_density_image(
        phases,
        moments_image,
        n_pixels,
        n_moments,
        spectral_image
    );
}


void bounded_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image)
{
    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image(
        normalized_moments_image,
        n_pixels,
        n_moments,
        mins, maxs,
        compressed_moments_image
    );

    bounded_decompress_moments_image(
        compressed_moments_image,
        n_pixels, n_moments,
        moments_image
    );

    compute_density_image(
        phases,
        moments_image,
        n_pixels,
        n_moments,
        spectral_image
    );
}


void unbounded_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image)
{
    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image(
        normalized_moments_image,
        n_pixels,
        n_moments,
        mins, maxs,
        compressed_moments_image
    );

    unbounded_decompress_moments_image(
        compressed_moments_image,
        n_pixels, n_moments,
        moments_image
    );

    compute_density_image(
        phases,
        moments_image,
        n_pixels,
        n_moments,
        spectral_image
    );
}


void unbounded_to_bounded_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image)
{
    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image(
        normalized_moments_image,
        n_pixels,
        n_moments,
        mins, maxs,
        compressed_moments_image
    );

    unbounded_to_bounded_decompress_moments_image(
        compressed_moments_image,
        n_pixels, n_moments,
        moments_image
    );

    compute_density_image(
        phases,
        moments_image,
        n_pixels,
        n_moments,
        spectral_image
    );
}

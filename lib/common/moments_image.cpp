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

#include "moments_image.h"
#include "moments.h"
#include "Util.h"

#include <cstring>
#include <Eigen/Core>


#define MIN_MAX

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
    assert(spectral_image.size() == phases.size() * n_pixels);

    moments_image.resize(n_pixels * n_moments);

    compute_moments_image(
        phases.data(),
        phases.size(),
        spectral_image.data(),
        n_pixels,
        n_moments,
        moments_image.data());

    assert(moments_image.size() == n_pixels * n_moments);
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
    assert(moments_image.size() == n_pixels * n_moments);

    compressed_moments_image.resize(moments_image.size());

    unbounded_compress_moments_image(
        moments_image.data(),
        n_pixels,
        n_moments,
        compressed_moments_image.data()
    );

    assert(compressed_moments_image.size() == n_pixels * n_moments);
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
    assert(compressed_moments_image.size() == n_pixels * n_moments);

    moments_image.resize(compressed_moments_image.size());

    unbounded_decompress_moments_image(
        compressed_moments_image.data(),
        n_pixels,
        n_moments,
        moments_image.data()
    );

    assert(moments_image.size() == n_pixels * n_moments);
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
    assert(moments_image.size() == n_pixels * n_moments);

    density_image.resize(phases.size() * n_pixels);

    compute_density_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        n_pixels,
        n_moments,
        density_image.data()
    );

    assert(density_image.size() == phases.size() * n_pixels);
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
    assert(moments_image.size() == n_pixels * n_moments);

    density_image.resize(phases.size() * n_pixels);

    bounded_compute_density_lagrange_image(
        phases.data(),
        phases.size(),
        moments_image.data(),
        n_pixels,
        n_moments,
        density_image.data()
    );

    assert(density_image.size() == phases.size() * n_pixels);
}


extern "C"
void normalize_moment_image_min_max(
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
            assert(min < max);

            dest[px * n_moments + m] = (src[px * n_moments + m] - min) / (max - min);
            assert(dest[px * n_moments + m] >= 0.);
            assert(dest[px * n_moments + m] <= 1.);
        }
    }
}


void normalize_moment_image_min_max(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    std::vector<double>& dest,
    std::vector<double>& mins,
    std::vector<double>& maxs)
{
    assert(src.size() == n_pixels * n_moments);

    dest.resize(src.size());
    mins.resize(n_moments - 1);
    maxs.resize(n_moments - 1);

    normalize_moment_image_min_max(
        src.data(),
        n_pixels, n_moments,
        dest.data(),
        mins.data(),
        maxs.data()
    );

    assert(dest.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
}

extern "C"
void denormalize_moment_image_min_max(
    const double src[],
    size_t n_pixels, size_t n_moments,
    const double mins[],
    const double maxs[],
    double dest[])
{
    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // DC component
        assert(!std::isnan(src[px * n_moments]));
        assert(!std::isinf(src[px * n_moments]));
        dest[px * n_moments] = src[px * n_moments];

        // AC components
        for (size_t m = 1; m < n_moments; m++) {
            const double& min = mins[m - 1];
            const double& max = maxs[m - 1];
            assert(min < max);

            const double normalized_moment = src[px * n_moments + m];

            assert(!std::isnan(normalized_moment));
            assert(!std::isinf(normalized_moment));
            dest[px * n_moments + m] = (max - min) * normalized_moment + min;
        }
    }
}

void denormalize_moment_image_min_max(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    std::vector<double>& dest)
{
    assert(src.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    dest.resize(src.size());

    denormalize_moment_image_min_max(
        src.data(),
        n_pixels, n_moments,
        mins.data(),
        maxs.data(),
        dest.data()
    );

    assert(dest.size() == n_pixels * n_moments);
}


extern "C"
void normalize_moment_image_stddev(
    const double src[],
    size_t n_pixels, size_t n_moments,
    double normalized[],
    double means[],
    double stddevs[])
{
    // Compute stats
    for (size_t m = 1; m < n_moments; m++) {
        double b[3] = {0};

        b[0] = n_pixels;

        for (size_t px = 0; px < n_pixels; px++) {
            const double px_value_1 = src[px * n_moments + m];
            const double px_value_2 = px_value_1 * px_value_1;

            b[1] += px_value_1;
            b[2] += px_value_2;
        }

        means[m - 1] = b[1] / b[0];
        stddevs[m - 1] = std::sqrt(b[0] * b[2] - b[1] * b[1]) / b[0];
    }

    // Rescale moments
    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // Ignore for 0th moment
        normalized[px * n_moments + 0] = src[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            normalized[px * n_moments + m] = (src[px * n_moments + m] - means[m - 1]) / (10. * stddevs[m - 1]);
            // assert(normalized[px * n_moments + m] >= -1.0);
            // assert(normalized[px * n_moments + m] <= 1.0);
            normalized[px * n_moments + m] = Util::clamp(0.5 * normalized[px * n_moments + m] + 0.5, 0.01, 0.99);
            // normalized[px * n_moments + m] = Util::clamp(normalized[px * n_moments + m], -0.99, 0.99);
        }
    }
}

void normalize_moment_image_stddev(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    std::vector<double>& dest,
    std::vector<double>& means,
    std::vector<double>& stddevs)
{
    assert(src.size() == n_pixels * n_moments);

    dest.resize(src.size());
    means.resize(n_moments - 1);
    stddevs.resize(n_moments - 1);

    normalize_moment_image_stddev(
        src.data(),
        n_pixels, n_moments,
        dest.data(),
        means.data(),
        stddevs.data()
    );

    assert(dest.size() == n_pixels * n_moments);
    assert(means.size() == n_moments - 1);
    assert(stddevs.size() == n_moments - 1);
}


extern "C"
void denormalize_moment_image_stddev(
    const double normalized[],
    size_t n_pixels, size_t n_moments,
    const double means[],
    const double stddevs[],
    double dest[])
{
    #pragma omp parallel for
    for (size_t px = 0; px < n_pixels; px++) {
        // Ignore for 0th moment
        dest[px * n_moments + 0] = normalized[px * n_moments + 0];

        for (size_t m = 1; m < n_moments; m++) {
            const double n = 2. * (normalized[px * n_moments + m] - 0.5);
            // const double n = normalized[px * n_moments + m];
            dest[px * n_moments + m] = n * 10. * stddevs[m - 1] + means[m - 1];
        }
    }
}

void denormalize_moment_image_stddev(
    const std::vector<double>& src,
    size_t n_pixels, size_t n_moments,
    const std::vector<double>& means,
    const std::vector<double>& stddevs,
    std::vector<double>& dest)
{
    assert(src.size() == n_pixels * n_moments);
    assert(means.size() == n_moments - 1);
    assert(stddevs.size() == n_moments - 1);

    dest.resize(src.size());

    denormalize_moment_image_stddev(
        src.data(),
        n_pixels, n_moments,
        means.data(),
        stddevs.data(),
        dest.data()
    );

    assert(dest.size() == n_pixels * n_moments);
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
    std::vector<double>& maxs,
    bool normalize_image)
{
    assert(spectral_image.size() == wavelengths.size() * n_pixels);

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

    if (normalize_image) {
        normalize_moment_image_min_max(
            moments_image,
            n_pixels, n_moments,
            normalized_moments_image,
            mins, maxs
        );
    } else {
        normalized_moments_image = moments_image;

        mins.resize(n_moments - 1);
        maxs.resize(n_moments - 1);

        for (size_t i = 0; i < n_moments - 1; i++) {
            mins[i] = 0.f;
            maxs[i] = 1.f;
        }
    }


    assert(normalized_moments_image.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
}


void linavg_compress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& spectral_image,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& normalized_moments_image,
    std::vector<double>& mins,
    std::vector<double>& maxs,
    bool normalize_image)
{
    assert(spectral_image.size() == wavelengths.size() * n_pixels);

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

    // Now divide every AC component by the average
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        for (size_t m = 1; m < n_moments; m++) {
            moments_image[i * n_moments + m] /= moments_image[i * n_moments];
        }
    }

    if (normalize_image) {
        normalize_moment_image_min_max(
            moments_image,
            n_pixels, n_moments,
            normalized_moments_image,
            mins, maxs
        );
    } else {
        normalized_moments_image = moments_image;

        mins.resize(n_moments - 1);
        maxs.resize(n_moments - 1);

        for (size_t i = 0; i < n_moments - 1; i++) {
            mins[i] = 0.f;
            maxs[i] = 1.f;
        }
    }

    assert(normalized_moments_image.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
}


void compress_spectral_image(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double> spectral_image,
    uint32_t width, uint32_t height,
    size_t n_moments,
    std::vector<double>& compressed_moments,
    std::vector<double>& mins, std::vector<double>& maxs,
    bool normalize_image,
    double& timing)
{
    const size_t n_pixels = width * height;

    auto clock_start = std::chrono::steady_clock::now();

    switch (method) {
        case LINEAR:
            linear_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs,
                normalize_image
            );
            break;

        case LINAVG:
            linavg_compress_spectral_image(
                wavelengths, spectral_image,
                n_pixels, n_moments,
                compressed_moments,
                mins, maxs,
                normalize_image
            );
            break;
    }

    auto clock_end = std::chrono::steady_clock::now();
    timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();

    assert(compressed_moments.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);
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
    assert(normalized_moments_image.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image_min_max(
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

    assert(spectral_image.size() == wavelengths.size() * n_pixels);
}


void linavg_decompress_spectral_image(
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image)
{
    assert(normalized_moments_image.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    std::vector<double> phases;
    std::vector<double> compressed_moments_image;
    std::vector<double> moments_image;

    wavelengths_to_phases(wavelengths, phases);

    denormalize_moment_image_min_max(
        normalized_moments_image,
        n_pixels,
        n_moments,
        mins, maxs,
        moments_image
    );

    // Multiply AC component by the average
    #pragma omp parallel for
    for (size_t i = 0; i < n_pixels; i++) {
        for (size_t m = 1; m < n_moments; m++) {
            moments_image[i * n_moments + m] *= moments_image[i * n_moments];
        }
    }

    compute_density_image(
        phases,
        moments_image,
        n_pixels,
        n_moments,
        spectral_image
    );

    assert(spectral_image.size() == wavelengths.size() * n_pixels);
}


void decompress_spectral_image(
    SpectralCompressionType method,
    const std::vector<double>& wavelengths,
    const std::vector<double>& normalized_moments_image,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    size_t n_pixels,
    size_t n_moments,
    std::vector<double>& spectral_image,
    double& timing)
{
    assert(normalized_moments_image.size() == n_pixels * n_moments);
    assert(mins.size() == n_moments - 1);
    assert(maxs.size() == n_moments - 1);

    auto clock_start = std::chrono::steady_clock::now();

    switch (method) {
        case LINEAR:
            linear_decompress_spectral_image(
                wavelengths,
                normalized_moments_image,
                mins, maxs,
                n_pixels, n_moments,
                spectral_image
            );
            break;

        case LINAVG:
            linavg_decompress_spectral_image(
                wavelengths,
                normalized_moments_image,
                mins, maxs,
                n_pixels, n_moments,
                spectral_image
            );
            break;
    }

    auto clock_end = std::chrono::steady_clock::now();
    timing = std::chrono::duration<double, std::milli>(clock_end - clock_start).count();

    assert(spectral_image.size() == wavelengths.size() * n_pixels);
}

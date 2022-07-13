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

#pragma once

#include <cstddef>
#include <vector>

void compute_moments_image(
    const std::vector<float>& phases,
    const std::vector<float>& spectral_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& moments_image
);


void compute_moments_image(
    const float phases[],
    size_t n_phases,
    const float spectral_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[]
);


void compress_moments_image(
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float compressed_moments_image[]
);


void compress_moments_image(
    const std::vector<float>& moments_image,
    size_t width, size_t height, 
    size_t n_moments,
    std::vector<float>& compressed_moments_image
);


void decompress_moments_image(
    const float compressed_moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float moments_image[]
);


void decompress_moments_image(
    const std::vector<float>& compressed_moments_image,
    size_t width, size_t height, 
    size_t n_moments,
    std::vector<float>& moments_image
);


void compute_density_image(
    const float phases[],
    size_t n_phases,
    const float moments_image[],
    size_t width, size_t height,
    size_t n_moments,
    float density_image[]
);


void compute_density_image(
    const std::vector<float>& phases,
    const std::vector<float>& moments_image,
    size_t width, size_t height,
    size_t n_moments,
    std::vector<float>& density_image
);
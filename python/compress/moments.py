#!/usr/bin/env python3

# Copyright (c) 2019, Christoph Peters
#               2022, Alban Fichet
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Karlsruhe Institute of Technology nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import numpy as np
import scipy

import compress.util as util


def wavelengths_to_phase(wavelengths: np.array, use_warp = False) -> np.array:
    if use_warp:
        warp_wavelengths = np.linspace(360.0, 830.0, 95)
        warp_phases = [
            -3.141592654, -3.130798150, -3.120409385, -3.109941320, -3.099293828, -3.088412485,
            -3.077521262, -3.065536918, -3.051316899, -3.034645062, -3.011566128, -2.977418908,
            -2.923394305, -2.836724273, -2.718861152, -2.578040247, -2.424210212, -2.263749529,
            -2.101393442, -1.939765501, -1.783840979, -1.640026237, -1.516709057, -1.412020736,
            -1.321787834, -1.241163202, -1.164633877, -1.087402539, -1.005125823, -0.913238812,
            -0.809194293, -0.691182886, -0.559164417, -0.416604255, -0.267950333, -0.117103760,
             0.033240425,  0.181302905,  0.326198172,  0.468115940,  0.607909250,  0.746661762,
             0.885780284,  1.026618635,  1.170211067,  1.316774725,  1.465751396,  1.615788104,
             1.764019821,  1.907870188,  2.044149116,  2.170497529,  2.285044356,  2.385874947,
             2.472481955,  2.546360090,  2.608833897,  2.660939270,  2.703842819,  2.739031973,
             2.767857527,  2.791481523,  2.810988348,  2.827540180,  2.842073425,  2.855038660,
             2.866280238,  2.876993433,  2.887283163,  2.897272230,  2.907046247,  2.916664580,
             2.926170139,  2.935595142,  2.944979249,  2.954494971,  2.964123227,  2.973828479,
             2.983585490,  2.993377472,  3.003193861,  3.013027335,  3.022872815,  3.032726596,
             3.042586022,  3.052449448,  3.062315519,  3.072183298,  3.082052181,  3.091921757,
             3.101791675,  3.111661747,  3.121531809,  3.131398522,  3.141592654
        ]

        return np.interp(wavelengths, warp_wavelengths, warp_phases)
    else:
        wl_min, wl_max = wavelengths[0], wavelengths[-1]
        phases = np.pi * (wavelengths - wl_min) / (wl_max - wl_min) - np.pi

        return phases


def signal_to_moments(phases: np.array, signal: np.array, n_moments: int) -> np.array:
    # We do not want to treat out of range values as zero but want to continue
    # the outermost value to the end of the integration range
    if phases[0] != -np.pi:
        phases = np.concatenate([[-np.pi], phases])
        signal = np.concatenate([[signal[0]], signal])
    if phases[-1] != 0:
        phases = np.concatenate([phases, [0]])
        signal = np.concatenate([signal, [signal[-1]]])

    trigonometric_moments = np.zeros((n_moments,), dtype=complex)

    phases_1 = phases[:-1]
    signal_1 = signal[:-1]

    phases_2 = phases[1:]
    signal_2 = signal[1:]

    d_signal = signal_2 - signal_1
    d_phases = phases_2 - phases_1

    a_k = d_signal / d_phases
    b_k = signal_1 - a_k * phases_1

    k = np.arange(0, n_moments)
    pp, kk = np.meshgrid(phases, k) # [m, p]
    exp = np.exp(-1j * pp * kk)

    cm_summand = (a_k / kk[1:, 1:] + 1j * b_k) / kk[1:, 1:]

    m2 = (cm_summand + a_k * 1j * phases_2 / kk[1:, 1:]) * exp[1:, 1:]
    m1 = (cm_summand + a_k * 1j * phases_1 / kk[1:, 1:]) * exp[1:, :-1]

    trigonometric_moments[1:] = np.sum(m2 - m1, axis=1)
    trigonometric_moments[0]  = np.sum(0.5 * a_k * phases_2 ** 2 + b_k * phases_2)
    trigonometric_moments[0] -= np.sum(0.5 * a_k * phases_1 ** 2 + b_k * phases_1)
    # trigonometric_moments[0] = np.sum(
    #     0.5 * a_k * (phases_2**2 - phases_1**2) + b_k * (phases_2 - phases_1)
    # )
    # trigonometric_moments[0] = np.sum((signal_2 - signal_1) * (phases_2 + phases_1) / 2 + b_k * (phases_2 - phases_1))

    # trigonometric_moments[0] = np.sum(
    #     d_signal * (phases_1 + phases_2) / 2
    #     + signal_1 * d_phases
    #     + d_signal * phases_1
    # )

    return trigonometric_moments.real / np.pi


def get_basis_signal_to_moments(phases: np.array) -> np.array:
    n_moments = len(phases)

    mat_signal_moments = np.zeros((n_moments, n_moments))

    signal = np.eye(n_moments)

    for i in range(n_moments):
        mat_signal_moments[i, :] = signal_to_moments(phases, signal[:, i], n_moments)

    return mat_signal_moments


def get_basis_moment_to_signal(phases: np.array) -> np.array:
    mat_moments_signal = get_basis_signal_to_moments(phases)

    return np.linalg.inv(mat_moments_signal)


def compute_density(phases: np.array, moments: np.array) -> np.array:
    n_moments = len(moments)

    e = np.zeros(len(moments))
    e[0] = 1

    evaluation_polynomial = scipy.linalg.solve_toeplitz(
        (moments / (2 * np.pi),
        np.conj(moments / (2 * np.pi))), e)

    k = np.arange(0, n_moments)
    pp, kk = np.meshgrid(phases, k)

    c = np.exp(-1j * pp * kk) / (2 * np.pi)

    num = evaluation_polynomial[0]
    denum = np.abs(evaluation_polynomial @ c)**2

    f = 1 / (2 * np.pi) * num / denum

    return f


def bounded_moments_to_exponential_moments(moments: np.array) -> np.array:
    exponential_moments = np.zeros_like(moments, dtype=complex)

    gamma_0_prime = 1/(4 * np.pi) * np.exp(1j * np.pi * (moments[0] - 1/2))
    exponential_moments[0] = 2 * np.real(gamma_0_prime)

    for l in range(1, exponential_moments.shape[0]):
        # sum = complex(0)
        # for j in range(1, l):
        #     sum += (l - j) * exponential_moments[j] * moments[l - j]
        sum = np.sum(np.arange(l - 1, 0, -1) * exponential_moments[1:l] * moments[l - 1:0:-1])

        exponential_moments[l] = (2 * np.pi * 1j / l) * (l * gamma_0_prime * moments[l] + sum)

    exponential_moments[0] = gamma_0_prime

    return exponential_moments


def bounded_exponential_moments_to_moments(exponential_moments: np.array) -> np.array:
    gamma_0_prime = exponential_moments[0]

    moments = np.zeros_like(exponential_moments, dtype=complex)
    moments[0] = 1/2 + np.log(4 * np.pi * gamma_0_prime) / (1j * np.pi)
    exponential_moments[0] = 2 * np.real(gamma_0_prime)

    for l in range(1, exponential_moments.shape[0]):
        # sum = complex(0)
        # for j in range(1, l):
        #     sum += (l - j) * exponential_moments[j] * moments[l - j]
        sum = np.sum(np.arange(l - 1, 0, -1) * exponential_moments[1:l] * moments[l - 1:0:-1])

        moments[l] = (
            exponential_moments[l] / (2j * np.pi * gamma_0_prime)
            - sum / (l * gamma_0_prime)
        )

    return np.real(moments)


def bounded_compress_real_trigonometric_moments(bounded_trigonometric_moments):
    """
    Maps the given real, bounded trigonometric moments to a representation that
    holds the same information using scalars in [-1,1], which are more suitable
    for compression through quantization.
    """
    exponential_moments = bounded_moments_to_exponential_moments(bounded_trigonometric_moments)

    toeplitz_first_column    = exponential_moments / (2 * np.pi)
    toeplitz_first_column[0] = 2.0 * toeplitz_first_column[0].real

    dots = util.get_levinson_dots(toeplitz_first_column)

    compressed = (dots * np.abs(exponential_moments[0]) / (1.0j * exponential_moments[0])).real
    compressed[0] = bounded_trigonometric_moments[0]

    return compressed


def unbounded_compress_real_trigonometric_moments(trigonometric_moments):
    """
    Like compress_real_bounded_trigonometric_moments() but for unbounded moments
    """
    toeplitz_first_column = trigonometric_moments
    dots = util.get_levinson_dots(toeplitz_first_column)
    compressed = dots
    compressed[0] = trigonometric_moments[0]
    return compressed


def unbounded_to_bounded_compress_real_trigonometric_moments(trigonometric_moments: np.array) -> np.array:
    bounded_trigonometric_moments = trigonometric_moments / (len(trigonometric_moments) * trigonometric_moments[0])

    exponential_moments = bounded_moments_to_exponential_moments(bounded_trigonometric_moments)

    toeplitz_first_column    = exponential_moments / (2 * np.pi)
    toeplitz_first_column[0] = 2.0 * toeplitz_first_column[0].real

    dots = util.get_levinson_dots(toeplitz_first_column)

    compressed = (dots * np.abs(exponential_moments[0]) / (1.0j * exponential_moments[0])).real
    compressed[0] = trigonometric_moments[0]

    return compressed


def bounded_decompress_real_trigonometric_moments(compressed):
    """Inverse mapping of bounded_compress_real_trigonometric_moments()."""
    exp_0 = 0.25 / np.pi * np.exp(np.pi * 1.0j * (compressed[0] - 0.5))

    dots = np.zeros_like(compressed, dtype=complex)
    dots[1:] = compressed[1:] * (1.0j * exp_0 / np.abs(exp_0))

    toeplitz_first_column = util.run_levinson_from_dots(np.real(exp_0) / np.pi, dots)

    exponential_moments = toeplitz_first_column * 2.0 * np.pi
    exponential_moments[0] = exp_0

    return np.real(bounded_exponential_moments_to_moments(exponential_moments))


def unbounded_decompress_real_trigonometric_moments(compressed):
    """Inverse mapping of unbounded_compress_real_trigonometric_moments()."""
    dots = np.zeros_like(compressed, dtype=complex)
    dots[1:] = compressed[1:]
    toeplitz_first_column = util.run_levinson_from_dots(compressed[0], dots)

    return toeplitz_first_column


def unbounded_to_bounded_decompress_real_trigonometric_moments(compressed_m):
    compressed = compressed_m
    compressed[0] = 1. / len(compressed)

    """Inverse mapping of bounded_compress_real_trigonometric_moments()."""
    exp_0 = 0.25 / np.pi * np.exp(np.pi * 1.0j * (compressed[0] - 0.5))

    dots = np.zeros_like(compressed, dtype=complex)
    dots[1:] = compressed[1:] * (1.0j * exp_0 / np.abs(exp_0))

    toeplitz_first_column = util.run_levinson_from_dots(np.real(exp_0) / np.pi, dots)

    exponential_moments = toeplitz_first_column * 2.0 * np.pi
    exponential_moments[0] = exp_0

    moments = np.real(bounded_exponential_moments_to_moments(exponential_moments))
    moments[0] = compressed_m[0]

    return moments


def bounded_forward(moment_image: np.array):
    w, h, n_moments = moment_image.shape

    moments = np.real(moment_image.reshape((w * h, n_moments)))

    compressed_moments  = np.zeros_like(moments)

    for i in range(moments.shape[0]):
        compressed_moments[i, :] = bounded_compress_real_trigonometric_moments(moments[i, :])

    normalized_moments, mins, maxs = util.normalize(compressed_moments)

    return normalized_moments.reshape((w, h, n_moments)), mins, maxs


def unbounded_forward(moment_image: np.array):
    w, h, n_moments = moment_image.shape

    moments = np.real(moment_image.reshape((w * h, n_moments)))

    compressed_moments  = np.zeros_like(moments)

    for i in range(moments.shape[0]):
        compressed_moments[i, :] = unbounded_compress_real_trigonometric_moments(moments[i, :])

    normalized_moments, mins, maxs = util.normalize(compressed_moments)

    return normalized_moments.reshape((w, h, n_moments)), mins, maxs


def unbounded_to_bounded_forward(moment_image: np.array):
    w, h, n_moments = moment_image.shape

    moments = np.real(moment_image.reshape((w * h, n_moments)))

    compressed_moments  = np.zeros_like(moments)

    for i in range(moments.shape[0]):
        m0 = moments[i, 0]
        moments[i, :] = moments[i, :] / (moments.shape[1] * m0)
        compressed_moments[i, :] = bounded_compress_real_trigonometric_moments(moments[i, :])
        compressed_moments[i, 0] = m0

    normalized_moments, mins, maxs = util.normalize(compressed_moments)

    return normalized_moments.reshape((w, h, n_moments)), mins, maxs


def bounded_backward(inv_base: np.array, normalized_moments_image: np.array, mins: np.array, maxs: np.array) -> np.array:
    w, h, n_moments = normalized_moments_image.shape

    normalized_moments = normalized_moments_image.reshape((w * h, n_moments))
    compressed_moments = util.denormalize(normalized_moments, mins, maxs)

    moments = np.zeros_like(compressed_moments)

    for i in range(moments.shape[0]):
        moments[i, :] = bounded_decompress_real_trigonometric_moments(compressed_moments[i, :])

    signals = moments @ inv_base

    return np.real(signals.reshape((w, h, n_moments)))


def unbounded_backward(inv_base: np.array, normalized_moments_image: np.array, mins: np.array, maxs: np.array) -> np.array:
    w, h, n_moments = normalized_moments_image.shape

    normalized_moments = normalized_moments_image.reshape((w * h, n_moments))
    compressed_moments = util.denormalize(normalized_moments, mins, maxs)

    moments = np.zeros_like(compressed_moments)

    for i in range(moments.shape[0]):
        moments[i, :] = unbounded_decompress_real_trigonometric_moments(compressed_moments[i, :])

    signals = moments @ inv_base

    return np.real(signals.reshape((w, h, n_moments)))


def unbounded_to_bounded_backward(inv_base: np.array, normalized_moments_image: np.array, mins: np.array, maxs: np.array) -> np.array:
    w, h, n_moments = normalized_moments_image.shape

    normalized_moments = normalized_moments_image.reshape((w * h, n_moments))
    compressed_moments = util.denormalize(normalized_moments, mins, maxs)

    moments = np.zeros_like(compressed_moments)

    for i in range(moments.shape[0]):
        m0 = compressed_moments[i, 0]
        compressed_moments[i, 0] = 1/moments.shape[1]
        moments[i, :] = bounded_decompress_real_trigonometric_moments(compressed_moments[i, :])
        moments[i, 1:] *= moments.shape[1] * m0
        moments[i, 0] = m0

    signals = moments @ inv_base

    return np.real(signals.reshape((w, h, n_moments)))
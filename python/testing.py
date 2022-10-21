#!/usr/bin/env python3

import numpy as np
import compress.moments as mom


def wavelength_to_phase(wavelength, wavelength_min, wavelength_max, mirror_signal, use_warp):
    """This function implements the various methods to map a wavelength in 
       nanometers to a phase between -pi and pi or -pi and 0 discussed in 
       Section 4.1 of the paper.
      \param wavelength The wavelength of light in nanometers (e.g. 360 to 830 
             for visible light). May be an array.
      \param wavelength_min,wavelength_max The minimal and maximal wavelengths 
             that may be passed through wavelength. Irrelevant if a warp is 
             used.
      \param mirror_signal Pass True to obtain phases between -pi and 0 only. 
             This also implies that a different warp is used.
      \param use_warp Pass False to use a linear mapping and True to use a 
             piecewise linear mapping that is optimized to make Fourier 
             coefficients resemble XYZ coefficients.
      \return An array of same shape as wavelength with entries between -pi and 
              pi if mirror_signal is False or -pi and 0 otherwise."""
    if use_warp:
        warp_wavelengths = np.linspace(360.0, 830.0, 95)
        if mirror_signal:
            warp_phases = [-3.141592654, -3.141592654, -3.141592654, -3.141592654, -3.141591857, -3.141590597,
                           -3.141590237, -3.141432053, -3.140119041, -3.137863071, -3.133438967, -3.123406739,
                           -3.106095749, -3.073470612, -3.024748900, -2.963566246, -2.894461907, -2.819659701,
                           -2.741784136, -2.660533432, -2.576526605, -2.490368187, -2.407962868, -2.334138406,
                           -2.269339880, -2.213127747, -2.162806279, -2.114787412, -2.065873394, -2.012511127,
                           -1.952877310, -1.886377224, -1.813129945, -1.735366957, -1.655108108, -1.573400329,
                           -1.490781436, -1.407519056, -1.323814008, -1.239721795, -1.155352390, -1.071041833,
                           -0.986956525, -0.903007113, -0.819061538, -0.735505101, -0.653346027, -0.573896987,
                           -0.498725202, -0.428534515, -0.363884284, -0.304967687, -0.251925536, -0.205301867,
                           -0.165356255, -0.131442191, -0.102998719, -0.079687644, -0.061092401, -0.046554594,
                           -0.035419229, -0.027113640, -0.021085743, -0.016716885, -0.013468661, -0.011125245,
                           -0.009497032, -0.008356318, -0.007571826, -0.006902676, -0.006366945, -0.005918355,
                           -0.005533442, -0.005193920, -0.004886397, -0.004601975, -0.004334090, -0.004077698,
                           -0.003829183, -0.003585923, -0.003346286, -0.003109231, -0.002873996, -0.002640047,
                           -0.002406990, -0.002174598, -0.001942639, -0.001711031, -0.001479624, -0.001248405,
                           -0.001017282, -0.000786134, -0.000557770, -0.000332262, 0.000000000]
        else:
            warp_phases = [-3.141592654, -3.130798150, -3.120409385, -3.109941320, -3.099293828, -3.088412485,
                           -3.077521262, -3.065536918, -3.051316899, -3.034645062, -3.011566128, -2.977418908,
                           -2.923394305, -2.836724273, -2.718861152, -2.578040247, -2.424210212, -2.263749529,
                           -2.101393442, -1.939765501, -1.783840979, -1.640026237, -1.516709057, -1.412020736,
                           -1.321787834, -1.241163202, -1.164633877, -1.087402539, -1.005125823, -0.913238812,
                           -0.809194293, -0.691182886, -0.559164417, -0.416604255, -0.267950333, -0.117103760,
                           0.033240425, 0.181302905, 0.326198172, 0.468115940, 0.607909250, 0.746661762, 0.885780284,
                           1.026618635, 1.170211067, 1.316774725, 1.465751396, 1.615788104, 1.764019821, 1.907870188,
                           2.044149116, 2.170497529, 2.285044356, 2.385874947, 2.472481955, 2.546360090, 2.608833897,
                           2.660939270, 2.703842819, 2.739031973, 2.767857527, 2.791481523, 2.810988348, 2.827540180,
                           2.842073425, 2.855038660, 2.866280238, 2.876993433, 2.887283163, 2.897272230, 2.907046247,
                           2.916664580, 2.926170139, 2.935595142, 2.944979249, 2.954494971, 2.964123227, 2.973828479,
                           2.983585490, 2.993377472, 3.003193861, 3.013027335, 3.022872815, 3.032726596, 3.042586022,
                           3.052449448, 3.062315519, 3.072183298, 3.082052181, 3.091921757, 3.101791675, 3.111661747,
                           3.121531809, 3.131398522, 3.141592654]
        return np.interp(wavelength, warp_wavelengths, warp_phases)
    else:
        normalized_wavelengths = (wavelength - wavelength_min) / (wavelength_max - wavelength_min)
        return ((1.0 if mirror_signal else 2.0) * np.pi) * normalized_wavelengths - np.pi

def compute_trigonometric_moments(phase, signal, moment_count, mirror_signal):
    """This function computes trigonometric moments for the given 2 pi periodic 
       signal. If the signal is bounded between zero and one, they will be 
       bounded trigonometric moments. phase has to be a sorted array, signal 
       must have same shape to provide corresponding values. The function uses 
       linear interpolation between the sample points. The values of the first 
       and last sample point are repeated towards the end of the domain. If 
       mirror_signal is True, the signal is mirrored at zero, which means that 
       the moments are real.
      \return An array of shape (moment_count+1,) where entry j is the j-th 
              trigonometric moment."""
    # We do not want to treat out of range values as zero but want to continue 
    # the outermost value to the end of the integration range
    phase = np.concatenate([[-np.pi], phase, [0.0 if mirror_signal else np.pi]])
    signal = np.concatenate([[signal[0]], signal, [signal[-1]]])
    # Now we are ready to integrate
    sample_count = phase.size
    trigonometric_moments = np.zeros((moment_count + 1,), dtype=np.complex)
    for l in range(sample_count - 1):
        if phase[l] >= phase[l + 1]:
            continue
        gradient = (signal[l + 1] - signal[l]) / (phase[l + 1] - phase[l])
        y_intercept = signal[l] - gradient * phase[l]
        for j in range(1, moment_count + 1):
            common_summands = gradient * 1.0 / j ** 2 + y_intercept * 1.0j / j
            trigonometric_moments[j] += (common_summands + gradient * 1.0j * j * phase[l + 1] / j ** 2) * np.exp(
                -1.0j * j * phase[l + 1])
            trigonometric_moments[j] -= (common_summands + gradient * 1.0j * j * phase[l] / j ** 2) * np.exp(
                -1.0j * j * phase[l])
        trigonometric_moments[0] += 0.5 * gradient * phase[l + 1] ** 2 + y_intercept * phase[l + 1]
        trigonometric_moments[0] -= 0.5 * gradient * phase[l] ** 2 + y_intercept * phase[l]
    trigonometric_moments *= 0.5 / np.pi
    if mirror_signal:
        return 2.0 * trigonometric_moments.real
    else:
        return trigonometric_moments

# def compute_trigonometric_moments(phases, signal, moment_count, mirror_signal):
#     """This function computes trigonometric moments for the given 2 pi periodic 
#        signal. If the signal is bounded between zero and one, they will be 
#        bounded trigonometric moments. phase has to be a sorted array, signal 
#        must have same shape to provide corresponding values. The function uses 
#        linear interpolation between the sample points. The values of the first 
#        and last sample point are repeated towards the end of the domain. If 
#        mirror_signal is True, the signal is mirrored at zero, which means that 
#        the moments are real.
#       \return An array of shape (moment_count+1,) where entry j is the j-th 
#               trigonometric moment."""
#     # We do not want to treat out of range values as zero but want to continue 
#     # the outermost value to the end of the integration range
#     if phases[0] != -np.pi:        
#         phases = np.concatenate([[-np.pi], phases])
#         signal = np.concatenate([[signal[0]], signal])
#     if mirror_signal and phases[-1] != 0:
#         phases = np.concatenate([phases, [0]])
#         signal = np.concatenate([signal, [signal[-1]]])
#     elif not mirror_signal and phases[-1] != np.pi:
#         phases = np.concatenate([phases, [np.pi]])
#         signal = np.concatenate([signal, [signal[-1]]])

#     trigonometric_moments = np.zeros((moment_count + 1,), dtype=complex)

#     phases_1 = phases[:-1]
#     phases_2 = phases[1:]

#     signal_1 = signal[:-1]
#     signal_2 = signal[1:]

#     a_k = (signal_2 - signal_1) / (phases_2 - phases_1)
#     b_k = signal_1 - a_k * phases_1


#     k = np.arange(0, moment_count + 1)
#     pp, kk = np.meshgrid(phases, k) # [m, p]
#     exp = np.exp(-1j * pp * kk)

#     cm_summand = (a_k / kk[1:, 1:] + 1j * b_k) / kk[1:, 1:]

#     m2 = (cm_summand + a_k * 1j * phases_2 / kk[1:, 1:]) * exp[1:, 1:]
#     m1 = (cm_summand + a_k * 1j * phases_1 / kk[1:, 1:]) * exp[1:, :-1]

#     trigonometric_moments[1:] = np.sum(m2 - m1, axis=1)
#     trigonometric_moments[0] = np.sum((phases_2 - phases_1) * ((signal_2 - signal_1) / 2 + b_k))

#     if mirror_signal:
#         return trigonometric_moments.real / np.pi
#     else:
#         return trigonometric_moments / (2 * np.pi)

def bounded_trigonometric_moments_to_exponential_moments(bounded_trigonometric_moments):
    """This function applies the recurrence specified Proposition 2 to turn 
       bounded trigonometric moments into exponential moments.
      \return An array of same shape as bounded_trigonometric_moments. Entry 0 
              does not hold the zeroth exponential moment but its complex 
              counterpart. Take its real part times two to get the actual 
              zeroth exponential moment."""
    moment_count = bounded_trigonometric_moments.size - 1
    exponential_moments = np.zeros((moment_count + 1,), dtype=complex)
    exponential_moments[0] = 0.25 / np.pi * np.exp(np.pi * 1.0j * (bounded_trigonometric_moments[0] - 0.5))
    for l in range(1, moment_count + 1):
        for j in range(l):
            exponential_moments[l] += (l - j) * exponential_moments[j] * bounded_trigonometric_moments[l - j]
        exponential_moments[l] *= 2.0j * np.pi / l
    return exponential_moments


def exponential_moments_to_bounded_trigonometric_moments(exponential_moments):
    """
    Inverse of bounded_trigonometric_moments_to_exponential_moments().
    """
    bounded = np.zeros_like(exponential_moments)
    bounded[0] = 1.0 / np.pi * np.angle(exponential_moments[0]) + 0.5
    exp_0 = 0.25 / np.pi * np.exp(np.pi * 1.0j * (bounded[0] - 0.5))
    for l in range(1, exponential_moments.shape[0]):
        bounded[l] = 1.0 / (2.0j * np.pi * exp_0) * exponential_moments[l] - 1.0 / (l * exp_0) * np.sum(np.arange(l - 1, 0, -1) * exponential_moments[1:l] * bounded[l - 1:0:-1])
    return bounded


def get_levinson_dots(first_column):
    """
    :return: A vector holding the dot products that arise in each iteration of
        Levinson's algorithm.
    """
    solution = np.zeros_like(first_column)
    solution[0] = 1.0 / first_column[0]
    dots = np.zeros_like(first_column)
    for j in range(1, first_column.shape[0]):
        dots[j] = np.dot(solution[0:j], first_column[j:0:-1])
        solution[0:j + 1] = (solution[0:j + 1] - dots[j] * np.conj(solution[j::-1])) / (1 - np.abs(dots[j]) ** 2)
    return dots


def run_levinson_from_dots(diagonal_entry, dots):
    """
    Runs Levinson's algorithm using the given intermediate dot products and the
    diagonal entry. Outputs the corresponding input vector and the result of
    Levinson's algorithm.
    """
    solution = np.zeros_like(dots)
    solution[0] = 1.0 / diagonal_entry
    first_column = np.zeros_like(dots)
    first_column[0] = diagonal_entry
    for j in range(1, dots.shape[0]):
        radius = 1.0 / solution[0].real
        center = -radius * np.dot(solution[1:j], first_column[j - 1:0:-1])
        first_column[j] = center + radius * dots[j]
        solution[0:j + 1] = (solution[0:j + 1] - dots[j] * np.conj(solution[j::-1])) / (1 - np.abs(dots[j]) ** 2)
    return first_column, solution


def compress_real_bounded_trigonometric_moments(bounded_trigonometric_moments):
    """
    Maps the given real, bounded trigonometric moments to a representation that
    holds the same information using scalars in [-1,1], which are more suitable
    for compression through quantization.
    """
    exponential_moments = bounded_trigonometric_moments_to_exponential_moments(bounded_trigonometric_moments)
    toeplitz_first_column = exponential_moments * 0.5 / np.pi
    toeplitz_first_column[0] = 2.0 * toeplitz_first_column[0].real
    dots = get_levinson_dots(toeplitz_first_column)
    compressed = (dots * np.abs(exponential_moments[0]) / (1.0j * exponential_moments[0])).real
    compressed[0] = bounded_trigonometric_moments[0]
    return compressed


def decompress_real_bounded_trigonometric_moments(compressed):
    """Inverse mapping of compress_real_bounded_trigonometric_moments()."""
    exp_0 = 0.25 / np.pi * np.exp(np.pi * 1.0j * (compressed[0] - 0.5))
    dots = np.zeros_like(compressed, dtype=complex)
    dots[1:] = compressed[1:] * (1.0j * exp_0 / np.abs(exp_0))
    toeplitz_first_column, evaluation_polynomial = run_levinson_from_dots(exp_0.real / np.pi, dots)
    exponential_moments = toeplitz_first_column * 2.0 * np.pi
    exponential_moments[0] = exp_0
    return exponential_moments_to_bounded_trigonometric_moments(exponential_moments).real


def compress_real_trigonometric_moments(trigonometric_moments):
    """
    Like compress_real_bounded_trigonometric_moments() but for unbounded moments
    """
    toeplitz_first_column = trigonometric_moments
    dots = get_levinson_dots(toeplitz_first_column)
    compressed = dots
    compressed[0] = trigonometric_moments[0]
    return compressed


def decompress_real_trigonometric_moments(compressed):
    """Inverse mapping of compress_real_trigonometric_moments()."""
    dots = np.zeros_like(compressed, dtype=complex)
    dots[1:] = compressed[1:]
    toeplitz_first_column, evaluation_polynomial = run_levinson_from_dots(compressed[0], dots)
    return toeplitz_first_column


def get_error(arr1, arr2):
    res = np.sum((arr1 - arr2)**2)

    return res


import struct
import matplotlib.pyplot as plt

def load_array(filename: str) -> np.array:
    with open(filename, 'rb') as f_c:
        data = f_c.read()
        sz, = struct.unpack_from('I', data)

        c_data = struct.unpack_from(
            str(sz) + 'd',
            data,
            offset=struct.calcsize('I'))

    # print(c_data)
    return np.asarray(c_data)


def main():
    wavelengths = np.linspace(312.5, 811.5, 999)
    emission    = np.asarray(
        [0.03183, 0.02949, 0.02075, 0.00807, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00281,
         0.00439, 0.00257, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00005, 0.00032, 0.00027, 0.00060, 0.00142, 0.00190, 0.00558, 0.00773, 0.00783, 0.00437, 0.00076, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00047, 0.00176, 0.00185, 0.00095, 0.00000, 0.00000,
         0.00266, 0.00169, 0.00037, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00090, 0.00475, 0.00823, 0.00661, 0.00546, 0.00271, 0.00265, 0.00472, 0.00521, 0.00544, 0.00392, 0.00160,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00536, 0.02746, 0.06113, 0.09360, 0.11145, 0.10533, 0.07810,
         0.05022, 0.02810, 0.01515, 0.00646, 0.00001, 0.00000, 0.00000, 0.00017, 0.00205, 0.00191, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00111, 0.00480, 0.00682, 0.00464, 0.00110, 0.00000, 0.00002, 0.00029, 0.00020,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00024, 0.00000,
         0.00022, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002,
         0.00118, 0.00329, 0.00520, 0.00666, 0.00434, 0.00184, 0.00179, 0.00059, 0.00115, 0.00195, 0.00340, 0.00471,
         0.01317, 0.04988, 0.12346, 0.20524, 0.24420, 0.22047, 0.14463, 0.07435, 0.03571, 0.02632, 0.02591, 0.02307,
         0.01642, 0.00911, 0.00342, 0.00165, 0.00169, 0.00268, 0.00328, 0.00329, 0.00192, 0.00317, 0.00422, 0.00453,
         0.00546, 0.00413, 0.00508, 0.00678, 0.00753, 0.00674, 0.00624, 0.00508, 0.00509, 0.00593, 0.00637, 0.00656,
         0.00608, 0.00531, 0.00467, 0.00528, 0.00658, 0.00730, 0.00797, 0.00652, 0.00652, 0.00885, 0.01169, 0.01330,
         0.01241, 0.00971, 0.00920, 0.01072, 0.01075, 0.01258, 0.01346, 0.01516, 0.01670, 0.01500, 0.01622, 0.02043,
         0.02490, 0.02979, 0.05001, 0.13486, 0.36258, 0.66126, 0.87177, 0.83723, 0.60594, 0.32710, 0.13781, 0.06077,
         0.03594, 0.02657, 0.02265, 0.01907, 0.01775, 0.01916, 0.02089, 0.02335, 0.02384, 0.02206, 0.01974, 0.01775,
         0.01984, 0.02336, 0.02659, 0.02621, 0.02557, 0.02689, 0.02814, 0.02638, 0.02050, 0.01707, 0.01526, 0.01758,
         0.01839, 0.01852, 0.01807, 0.01876, 0.02014, 0.02036, 0.02159, 0.01983, 0.01828, 0.01855, 0.01811, 0.01768,
         0.01646, 0.01450, 0.01543, 0.01747, 0.01784, 0.01792, 0.01730, 0.01729, 0.01780, 0.01835, 0.01960, 0.01813,
         0.01523, 0.01225, 0.01116, 0.01358, 0.01481, 0.01464, 0.01427, 0.01450, 0.01595, 0.01605, 0.01544, 0.01418,
         0.01303, 0.01319, 0.01296, 0.01464, 0.01539, 0.01576, 0.01787, 0.02032, 0.01991, 0.01901, 0.01778, 0.01765,
         0.01899, 0.01974, 0.01968, 0.02081, 0.02057, 0.01812, 0.01486, 0.01408, 0.01631, 0.01952, 0.02394, 0.02725,
         0.03111, 0.03460, 0.03862, 0.04397, 0.05210, 0.06303, 0.07581, 0.09005, 0.10531, 0.12224, 0.13913, 0.15221,
         0.15977, 0.16024, 0.15716, 0.15020, 0.14103, 0.13206, 0.12378, 0.11772, 0.11313, 0.11063, 0.10707, 0.10022,
         0.09238, 0.08555, 0.08027, 0.07729, 0.07392, 0.07113, 0.06808, 0.06391, 0.06047, 0.05592, 0.05356, 0.04986,
         0.04574, 0.04123, 0.03665, 0.03348, 0.03017, 0.02693, 0.02394, 0.02124, 0.01992, 0.01856, 0.01863, 0.01795,
         0.01799, 0.01805, 0.01815, 0.01741, 0.01693, 0.01699, 0.01712, 0.01586, 0.01384, 0.01203, 0.01213, 0.01416,
         0.01523, 0.01541, 0.01526, 0.01498, 0.01456, 0.01330, 0.01246, 0.01241, 0.01309, 0.01443, 0.01411, 0.01512,
         0.01543, 0.01611, 0.01622, 0.01485, 0.01237, 0.01089, 0.01029, 0.01168, 0.01362, 0.01359, 0.01260, 0.01080,
         0.01076, 0.01122, 0.01149, 0.01112, 0.01209, 0.01326, 0.01394, 0.01389, 0.01420, 0.01471, 0.01493, 0.01390,
         0.01315, 0.01370, 0.01517, 0.01632, 0.01684, 0.01749, 0.02101, 0.02619, 0.03273, 0.03792, 0.04030, 0.04103,
         0.04012, 0.04190, 0.04794, 0.05849, 0.07639, 0.09780, 0.12232, 0.14726, 0.17818, 0.22662, 0.30577, 0.43080,
         0.59592, 0.76055, 0.87011, 0.88186, 0.82737, 0.76166, 0.72346, 0.74437, 0.85165, 1.06570, 1.28210, 1.32660,
         1.17590, 0.90671, 0.67591, 0.54440, 0.46602, 0.40838, 0.35627, 0.30960, 0.26661, 0.22940, 0.19965, 0.17487,
         0.15402, 0.13470, 0.11625, 0.09974, 0.08546, 0.07338, 0.06383, 0.05677, 0.05193, 0.04969, 0.04491, 0.03973,
         0.03454, 0.03274, 0.03187, 0.03158, 0.02883, 0.02584, 0.02381, 0.02323, 0.02310, 0.02224, 0.02066, 0.02031,
         0.02031, 0.02064, 0.01942, 0.01802, 0.01606, 0.01564, 0.01538, 0.01521, 0.01488, 0.01492, 0.01464, 0.01369,
         0.01246, 0.01223, 0.01135, 0.01102, 0.01003, 0.01086, 0.01199, 0.01359, 0.01600, 0.02732, 0.05904, 0.11278,
         0.16184, 0.18317, 0.17462, 0.17225, 0.19120, 0.21025, 0.20751, 0.18799, 0.16706, 0.15800, 0.16107, 0.17413,
         0.18906, 0.19824, 0.20043, 0.19760, 0.19154, 0.18733, 0.19436, 0.22298, 0.26762, 0.30149, 0.29477, 0.24851,
         0.18934, 0.14842, 0.12798, 0.12800, 0.13592, 0.15307, 0.17962, 0.21513, 0.25133, 0.27003, 0.25396, 0.21147,
         0.16101, 0.12577, 0.10811, 0.10260, 0.10423, 0.10871, 0.12116, 0.13824, 0.16119, 0.18220, 0.18892, 0.17842,
         0.15230, 0.12370, 0.10054, 0.08511, 0.07577, 0.07104, 0.07128, 0.07507, 0.07836, 0.08222, 0.08907, 0.10029,
         0.11541, 0.13450, 0.16321, 0.20835, 0.28525, 0.42854, 0.74318, 1.28260, 1.98160, 2.50960, 2.60340, 2.22100,
         1.68590, 1.25740, 1.00590, 0.86514, 0.75714, 0.64514, 0.52353, 0.41218, 0.32364, 0.26156, 0.21912, 0.19072,
         0.17270, 0.16387, 0.16293, 0.16424, 0.16667, 0.16742, 0.17001, 0.17361, 0.18034, 0.18947, 0.19862, 0.20592,
         0.21157, 0.21577, 0.22009, 0.22143, 0.21836, 0.21096, 0.20303, 0.19807, 0.19681, 0.20724, 0.23761, 0.28877,
         0.33451, 0.34214, 0.29759, 0.22037, 0.14536, 0.09589, 0.06834, 0.05638, 0.04923, 0.04675, 0.04636, 0.04507,
         0.04072, 0.03586, 0.03298, 0.03357, 0.03433, 0.03388, 0.03239, 0.02874, 0.02497, 0.02097, 0.02374, 0.02725,
         0.02870, 0.02767, 0.02611, 0.02723, 0.02679, 0.02596, 0.02709, 0.02802, 0.02793, 0.02496, 0.02355, 0.02643,
         0.03300, 0.04684, 0.06912, 0.09304, 0.10242, 0.09417, 0.07694, 0.06300, 0.05589, 0.05033, 0.04560, 0.03827,
         0.03365, 0.02854, 0.02634, 0.02628, 0.02984, 0.03226, 0.03216, 0.02992, 0.02733, 0.02439, 0.02408, 0.02443,
         0.02613, 0.03074, 0.04197, 0.05554, 0.06506, 0.05996, 0.04694, 0.03306, 0.02773, 0.02519, 0.02387, 0.02143,
         0.01990, 0.01292, 0.00916, 0.00956, 0.01480, 0.01783, 0.01691, 0.01702, 0.01404, 0.01028, 0.00776, 0.01151,
         0.01353, 0.01262, 0.00752, 0.00772, 0.00663, 0.00480, 0.00452, 0.00605, 0.00982, 0.00942, 0.00863, 0.00951,
         0.01099, 0.01282, 0.01298, 0.01573, 0.01831, 0.01877, 0.01796, 0.01318, 0.00856, 0.00709, 0.01268, 0.01888,
         0.02218, 0.02146, 0.02032, 0.02244, 0.03325, 0.04689, 0.05390, 0.05241, 0.04452, 0.03212, 0.01950, 0.01029,
         0.00849, 0.00685, 0.00676, 0.00815, 0.01425, 0.02679, 0.03881, 0.04250, 0.03391, 0.01795, 0.00639, 0.00254,
         0.00522, 0.00578, 0.00281, 0.00134, 0.00004, 0.00054, 0.00135, 0.00267, 0.00420, 0.00272, 0.00717, 0.01148,
         0.01327, 0.01113, 0.01344, 0.02051, 0.02905, 0.04216, 0.06547, 0.10035, 0.12991, 0.14735, 0.14673, 0.14957,
         0.15497, 0.15811, 0.13606, 0.10156, 0.07295, 0.07317, 0.09199, 0.11486, 0.12240, 0.10485, 0.07150, 0.04066,
         0.02214, 0.01333, 0.00693, 0.00168, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00083, 0.00621, 0.00822, 0.00297, 0.00053,
         0.00000, 0.00000, 0.00000, 0.00069, 0.00016, 0.00000, 0.00000, 0.00196, 0.00421, 0.00059, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00125, 0.00026,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00202, 0.02020, 0.02191, 0.01320, 0.00336, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00249, 0.00022, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00138, 0.00000, 0.00000, 0.00000, 0.00572, 0.00518, 0.00000, 0.00000,
         0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00119, 0.01094, 0.01152,
         0.01058, 0.02523, 0.04864])

    emission /= np.max(emission)
    print('signal bounds: [', np.min(emission), ',', np.max(emission), ']')

    n_moments = len(wavelengths)


    # Original implementation
    og_phases                = wavelength_to_phase(wavelengths, wavelengths.min(), wavelengths.max(), True, False)
    og_moments               = compute_trigonometric_moments(og_phases, emission, n_moments - 1, True)
    
    og_ub_compressed_moments = compress_real_trigonometric_moments(og_moments)
    og_d_ub_moments          = decompress_real_trigonometric_moments(og_ub_compressed_moments)

    og_b_compressed_moments  = compress_real_bounded_trigonometric_moments(og_moments)
    og_d_b_moments           = decompress_real_bounded_trigonometric_moments(og_b_compressed_moments)


    # Refactored implementation
    refactor_phases                = mom.wavelengths_to_phase(wavelengths)
    refactor_moments               = mom.signal_to_moments(refactor_phases, emission, n_moments)

    refactor_ub_compressed_moments = mom.unbounded_compress_real_trigonometric_moments(refactor_moments)
    refactor_d_ub_moments          = mom.unbounded_decompress_real_trigonometric_moments(refactor_ub_compressed_moments)

    refactor_b_compressed_moments  = mom.bounded_compress_real_trigonometric_moments(refactor_moments)
    refactor_d_b_moments           = mom.bounded_decompress_real_trigonometric_moments(refactor_b_compressed_moments)

    # b = mom.get_basis_moment_to_signal(refactor_phases)

    # refactor_signal = refactor_d_b_moments @ b

    print()
    print('------------------------------------')
    print('-- difference reference, refactor --')
    print('------------------------------------')
    print('              phases:', get_error(og_phases, refactor_phases))
    print('             moments:', get_error(og_moments, refactor_moments))
    print('  bounded compressed:', get_error(og_b_compressed_moments, refactor_b_compressed_moments))
    print('unbounded compressed:', get_error(og_ub_compressed_moments, refactor_ub_compressed_moments))
    print()
    print('------------------------------------')
    print('--   Comparison forward reverse   --')
    print('------------------------------------')
    print('  bounded original:', get_error(og_moments, og_d_b_moments))
    print('unbounded original:', get_error(og_moments, og_d_ub_moments))
    print()
    print('  bounded refactor:', get_error(refactor_moments, refactor_d_b_moments))
    print('unbounded refactor:', get_error(refactor_moments, refactor_d_ub_moments))
    print()
    print(np.min(og_b_compressed_moments), np.max(og_b_compressed_moments))
    print(np.min(og_ub_compressed_moments), np.max(og_ub_compressed_moments))
    print()
    print(np.min(refactor_b_compressed_moments), np.max(refactor_b_compressed_moments))
    print(np.min(refactor_ub_compressed_moments), np.max(refactor_ub_compressed_moments))

    print('-- C++ diff --')
    print(get_error(load_array('phases.dat'), og_phases))
    print(get_error(load_array('moments.dat'), og_moments))
    print(get_error(load_array('c_moments_b.dat'), refactor_b_compressed_moments))
    print(get_error(load_array('c_moments_ub.dat'), refactor_ub_compressed_moments))

if __name__ == '__main__':
    main()
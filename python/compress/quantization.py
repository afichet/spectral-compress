#!/usr/bin/env python3

import numpy as np
import moments


def quantize_dequantize(vec: np.array, index: int, n_bits: int):
    vec_q = vec.copy()
    
    q_value = (1 << n_bits) - 1
    vec_q[:, index] = np.round(vec[:, index] * q_value) / q_value

    return vec_q


def generate_bounded_quantization_curve(wavelengths: np.array, ref: np.array, n_bits: int, err_fun=rrmse):
    phases = moments.wavelength_to_phase(wavelengths)
    signal_moments, moments_signal = moments.build_base(phases)
    norm_moments, mins, maxs = moments.bounded_forward(signal_moments, ref)

    bits = np.zeros((norm_moments.shape[1],), dtype=np.uint8)

    bits[0] = 16
    bits[1] = n_bits

    # Determine error baseline
    q_norm_moments  = quantize_dequantize(norm_moments, 1, n_bits)
    backward_signal = moments.bounded_backward(moments_signal, q_norm_moments, mins, maxs)

    sz = backward_signal.shape[0] * backward_signal.shape[1]
    err_init = err_fun(wavelengths, ref, backward_signal)

    for i in range(2, norm_moments.shape[1]):
        bits[i] = bits[i - 1]
        
        for b in range(bits[i], 0, -1):
            # Decrease bitrate while keeping the error bellow our threshold
            q_norm_moments = quantize_dequantize(norm_moments, i, b)
            backward_signal = moments.bounded_backward(moments_signal, q_norm_moments, mins, maxs)

            err = np.linalg.norm(ref - backward_signal) / np.sqrt(sz)
            err = err_fun(wavelengths, ref, backward_signal)

            if err >= err_init:
                break
                
            bits[i] = b
    
    return bits
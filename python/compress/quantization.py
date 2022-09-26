#!/usr/bin/env python3

import numpy as np
import compress.moments as moments


def quantize_dequantize(vec: np.array, index: int, n_bits: int):
    vec_q = vec.copy()
    
    q_value = (1 << n_bits) - 1
    vec_q[:, index] = np.round(vec[:, index] * q_value) / q_value

    return vec_q


def rrmse(wl: np.array, s_ref: np.array, s_cmp: np.array) -> float:
    return np.sqrt(np.sum((s_ref - s_cmp)**2)) / np.average(s_ref)


def bounded_get_error_progressive_quantization(wavelengths: np.array, ref: np.array, n_bits: int, err_fun=rrmse) -> np.array:
    phases = moments.wavelengths_to_phase(wavelengths)
    
    signal_to_moments = moments.get_basis_signal_to_moments(phases)
    moments_to_signal = np.linalg.inv(signal_to_moments)

    norm_moments, mins, maxs = moments.bounded_forward(signal_to_moments, ref)

    err = np.zeros((ref.shape[1], ))

    for i in range(1, norm_moments.shape[1]):
        q_norm_moments = quantize_dequantize(norm_moments, i, n_bits)
        backward_signal = moments.bounded_backward(moments_to_signal, q_norm_moments, mins, maxs)

        err[i] = err_fun(wavelengths, ref, backward_signal)

    return err


def bounded_generate_quantization_curve(wavelengths: np.array, ref: np.array, n_bits: int, err_fun=rrmse) -> np.array:
    phases = moments.wavelengths_to_phase(wavelengths)

    signal_to_moments = moments.get_basis_signal_to_moments(phases)
    moments_to_signal = np.linalg.inv(signal_to_moments)

    norm_moments, mins, maxs = moments.bounded_forward(signal_to_moments, ref)

    bits = np.zeros((norm_moments.shape[1],), dtype=np.uint8)

    bits[0] = 16
    bits[1] = n_bits

    # Determine error baseline
    q_norm_moments  = quantize_dequantize(norm_moments, 1, n_bits)
    backward_signal = moments.bounded_backward(moments_to_signal, q_norm_moments, mins, maxs)

    sz = backward_signal.shape[0] * backward_signal.shape[1]
    err_init = err_fun(wavelengths, ref, backward_signal)

    for i in range(2, norm_moments.shape[1]):
        bits[i] = bits[i - 1]
        
        for b in range(bits[i], 0, -1):
            # Decrease bitrate while keeping the error bellow our threshold
            q_norm_moments = quantize_dequantize(norm_moments, i, b)
            backward_signal = moments.bounded_backward(moments_to_signal, q_norm_moments, mins, maxs)

            err = np.linalg.norm(ref - backward_signal) / np.sqrt(sz)
            err = err_fun(wavelengths, ref, backward_signal)

            if err >= err_init:
                break
                
            bits[i] = b
    
    return bits


def bounded_err_for_curve(phases: np.array, ref: np.array, curve: np.array, err_fun=rrmse):   
    signal_to_moments = moments.get_basis_signal_to_moments(phases)
    moments_to_signal = np.linalg.inv(signal_to_moments)

    norm_moments, mins, maxs = moments.bounded_forward(signal_to_moments, ref)

    quantized_norm_moments = norm_moments.copy()
    
    for i in range(1, norm_moments.shape[1]):
        quantized_norm_moments = quantize_dequantize(quantized_norm_moments, i, curve[i])

    backward_signal = moments.bounded_backward(moments_to_signal, quantized_norm_moments, mins, maxs)

    total_bits = np.sum(curve)
    err = err_fun(ref, backward_signal)

    return total_bits, err
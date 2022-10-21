#!/usr/bin/env python3

import unittest
import wave
import numpy as np
import compress.moments as mnt
from numpy.testing import assert_almost_equal

from test_data.sample_spectra import reflectance_wavelengths, reflectance, emission_wavelengths, emission


class TestMoments(unittest.TestCase):
    def test_signal_to_moment_to_signal(self):
        for wl, spectrum in [
            (reflectance_wavelengths, reflectance), 
            (emission_wavelengths, emission)]:
            phases = mnt.wavelengths_to_phase(wl)

            s_to_m = mnt.get_basis_signal_to_moments(phases)
            m_to_s = mnt.get_basis_moment_to_signal(phases)

            moments_fun = mnt.signal_to_moments(phases, spectrum, len(phases))
            moments_mat = spectrum @ s_to_m

            assert_almost_equal(moments_fun, moments_mat)

            signal = moments_fun @ m_to_s

            assert_almost_equal(signal, spectrum)


    def test_unbounded_compress(self):
        for wl, spectrum in [
            (reflectance_wavelengths, reflectance), 
            (emission_wavelengths, emission)]:
            phases = mnt.wavelengths_to_phase(wl)
            moments = mnt.signal_to_moments(phases, spectrum, len(phases))

            compressed_moments = mnt.unbounded_compress_real_trigonometric_moments(moments)
            decompressed = mnt.unbounded_decompress_real_trigonometric_moments(compressed_moments)

            assert_almost_equal(moments, decompressed)


    def test_bounded_compress(self):
        for wl, spectrum in [
            (reflectance_wavelengths, reflectance), 
            (emission_wavelengths, emission)]:
            if (np.max(spectrum) > 1):
                em = spectrum / np.max(spectrum)
            else:
                em = spectrum
            phases = mnt.wavelengths_to_phase(wl)
            moments = mnt.signal_to_moments(phases, em, len(phases))

            compressed_moments = mnt.bounded_compress_real_trigonometric_moments(moments)
            decompressed = mnt.bounded_decompress_real_trigonometric_moments(compressed_moments)

            assert_almost_equal(moments, decompressed)


    def test_unbounded_to_bounded_compress(self):
        for wl, spectrum in [
            (reflectance_wavelengths, reflectance), 
            (emission_wavelengths, emission)]:
            phases = mnt.wavelengths_to_phase(wl)
            moments = mnt.signal_to_moments(phases, spectrum, len(phases))

            compressed_moments = mnt.unbounded_to_bounded_compress_real_trigonometric_moments(moments)
            decompressed = mnt.unbounded_to_bounded_decompress_real_trigonometric_moments(compressed_moments)

            assert_almost_equal(moments, decompressed)


    # def test_full_pipeline(self):


if __name__ == '__main__':
    unittest.main()
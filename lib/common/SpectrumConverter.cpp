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


#include "SpectrumConverter.h"
#include "spectral_data.h"
#include "Util.h"

#include <cstring>
#include <cassert>
#include <cmath>
#include <functional>

SpectrumConverter::SpectrumConverter(bool emissiveSpectrum)
    : _emissiveSpectrum(emissiveSpectrum)
    , _cmfFirstWavelength_nm(CIE1931_2DEG_FIRST_WAVELENGTH_NM)
{
    if (!_emissiveSpectrum) {
        _illuminantFirstWavelenght_nm = D_65_FIRST_WAVELENGTH_NM;
        _illuminantSPD = std::vector<float>(std::begin(D_65_SPD), std::end(D_65_SPD));
    }

    _xyzCmfs[0] = std::vector<float>(
        std::begin(CIE1931_2DEG_X),
        std::end(CIE1931_2DEG_X));
    _xyzCmfs[1] = std::vector<float>(
        std::begin(CIE1931_2DEG_Y),
        std::end(CIE1931_2DEG_Y));
    _xyzCmfs[2] = std::vector<float>(
        std::begin(CIE1931_2DEG_Z),
        std::end(CIE1931_2DEG_Z));

    memcpy(&_xyzToRgb[0], XYZ_TO_SRGB_D65_MATRIX, 9 * sizeof(float));
}


SpectrumConverter::SpectrumConverter(
    const float &cmfFirstWavelength_nm,
    const std::array<std::vector<float>, 3> &xyzCmfs,
    const std::array<float, 9> xyzToRgb)
    : _cmfFirstWavelength_nm(cmfFirstWavelength_nm)
    , _xyzCmfs(xyzCmfs)
    , _xyzToRgb(xyzToRgb)
{}


size_t SpectrumConverter::cmfWavelengthIndex(float wavelength_nm) const
{
    assert(wavelength_nm >= firstWavelength());
    assert(wavelength_nm <= lastWavelength());

    float idx_f = float(wavelength_nm - _cmfFirstWavelength_nm);

    assert(idx_f >= 0);
    assert(idx_f < _xyzCmfs[0].size());

    return static_cast<size_t>(std::round(idx_f));
}


size_t SpectrumConverter::cmfWavelengthValue(size_t index) const
{
    // assert(index < _xyzCmfs[0].size());
    if (index >= _xyzCmfs[0].size()) {
        return 0;
    }

    return _cmfFirstWavelength_nm + index;
}


void SpectrumConverter::spectralImageToRGB(
    const std::vector<float>& wavelengths_nm,
    const std::vector<float>& spectral_image,
    size_t width, size_t height,
    std::vector<float>& rgb_image) const
{
    rgb_image.resize(3 * width * height);

    #pragma omp parallel for
    for (size_t i = 0; i < width * height; i++) {
        std::array<float, 3> rgb;

        spectrumToRGB(
            wavelengths_nm,
            &spectral_image[i * wavelengths_nm.size()],
            rgb);

        memcpy(&rgb_image[3 * i], &rgb[0], 3 * sizeof(float));
    }
}

void SpectrumConverter::spectrumToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *spectrum,
    std::array<float, 3> &XYZ) const
{
    if (_emissiveSpectrum) {
        emissiveSpectrumToXYZ(wavelengths_nm, spectrum, XYZ);
    } else {
        reflectiveSpectrumToXYZ(wavelengths_nm, spectrum, XYZ);
    }
}


void SpectrumConverter::spectrumToRGB(
    const std::vector<float> &wavelengths_nm,
    const float *spectrum,
    std::array<float, 3> &RGB) const
{
    std::array<float, 3> XYZ;
    spectrumToXYZ(wavelengths_nm, spectrum, XYZ);

    memset(&RGB[0], 0, 3 * sizeof(float));

    // Convert to RGB using the provided matrix
    for (size_t channel = 0; channel < 3; channel++) {
        for (size_t col = 0; col < 3; col++) {
            RGB[channel] += XYZ[col] * _xyzToRgb[3 * channel + col];
        }

        // Ensure RGB values are > 0
        RGB[channel] = std::max(RGB[channel], 0.F);
    }
}


void SpectrumConverter::spectraToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *reflectiveSpectrum,
    const float *emissiveSpectrum,
    std::array<float, 3> &XYZ) const
{
    std::array<float, 3> XYZ_refl;
    std::array<float, 3> XYZ_emissive;

    reflectiveSpectrumToXYZ(wavelengths_nm, reflectiveSpectrum, XYZ_refl);
    emissiveSpectrumToXYZ(wavelengths_nm, emissiveSpectrum, XYZ_emissive);

    for (size_t c = 0; c < 3; c++) {
        XYZ[c] = XYZ_refl[c] + XYZ_emissive[c];
    }
}


void SpectrumConverter::spectraToRGB(
    const std::vector<float> &wavelengths_nm,
    const float *reflectiveSpectrum,
    const float *emissiveSpectrum,
    std::array<float, 3> &RGB) const
{
    std::array<float, 3> XYZ;
    spectraToXYZ(wavelengths_nm, reflectiveSpectrum, emissiveSpectrum, XYZ);

    memset(&RGB[0], 0, 3 * sizeof(float));

    // Convert to RGB using the provided matrix
    for (size_t channel = 0; channel < 3; channel++) {
        for (size_t col = 0; col < 3; col++) {
            RGB[channel] += XYZ[col] * _xyzToRgb[3 * channel + col];
        }

        // Ensure RGB values are > 0
        RGB[channel] = std::max(RGB[channel], 0.F);
    }
}


void SpectrumConverter::spectrumToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *diagonal,
    const float *reradiation,
    std::array<float, 3> &XYZ) const
{
    memset(&XYZ[0], 0, 3 * sizeof(float));

    if (_emissiveSpectrum) {
        return spectrumToXYZ(wavelengths_nm, diagonal, XYZ);
    }

    if (wavelengths_nm.size() == 0) {
        return;
    }

    const float illuminant_last_wavelength = _illuminantFirstWavelenght_nm + _illuminantSPD.size() - 1;
    const float start_wavelength = std::max(
        std::max(_illuminantFirstWavelenght_nm, firstWavelength()),
        wavelengths_nm.front());
    const float end_wavelength = std::min(
        std::min(illuminant_last_wavelength, lastWavelength()),
        wavelengths_nm.back());

    // Early exit, selection out of range
    if (end_wavelength < start_wavelength) {
        return;
    }

    assert(start_wavelength <= end_wavelength);

    float normalisation_factor(0);

    std::function<float(size_t, size_t)> value = [diagonal, reradiation](size_t wi, size_t wo) {
        if (wi == wo) {
            return diagonal[wi];
        }

        return reradiation[Util::idxFromWavelengthIdx(wi, wo)];
    };

    for (size_t wl_idx_i = 0; wl_idx_i < wavelengths_nm.size() - 1;
         wl_idx_i++) {
        float wl_i_a = wavelengths_nm[wl_idx_i];
        float wl_i_b = wavelengths_nm[wl_idx_i + 1];

        // We have not reached yet the starting point
        if (start_wavelength > wl_i_b) {
            continue;
        }

        // We have finished the integration
        if (end_wavelength < wl_i_a) {
            break;
        }

        if (start_wavelength > wl_i_a) {
            wl_i_a = start_wavelength;
        }

        if (end_wavelength < wl_i_b) {
            wl_i_b = end_wavelength;
        }

        const size_t idx_illu_start = wl_i_a - _illuminantFirstWavelenght_nm;
        size_t idx_illu_end = wl_i_b - _illuminantFirstWavelenght_nm;

        // On last intervall we need to include the last
        // wavelength of the spectrum
        if (wl_idx_i == wavelengths_nm.size() - 2) {
            idx_illu_end = idx_illu_end + 1;
        }

        for (size_t wl_idx_o = wl_idx_i;
             wl_idx_o < wavelengths_nm.size() - 1;
             wl_idx_o++) {
            float wl_o_a = wavelengths_nm[wl_idx_o];
            float wl_o_b = wavelengths_nm[wl_idx_o + 1];

            // We have not reached yet the starting point
            if (start_wavelength > wl_o_b) {
                continue;
            }

            // We have finished the integration
            if (end_wavelength < wl_o_a) {
                break;
            }

            if (start_wavelength > wl_o_a) {
                wl_o_a = start_wavelength;
            }

            if (end_wavelength < wl_o_b) {
                wl_o_b = end_wavelength;
            }

            const size_t idx_cmf_start = cmfWavelengthIndex(wl_o_a);
            size_t idx_cmf_end = cmfWavelengthIndex(wl_o_b);

            // On last intervall we need to include the last
            // wavelength of the spectrum
            if (wl_idx_o == wavelengths_nm.size() - 2) {
                idx_cmf_end = idx_cmf_end + 1;
            }

            for (size_t idx_illu = idx_illu_start; idx_illu < idx_illu_end;
                 idx_illu++) {
                const float curr_wl_i = idx_illu + _illuminantFirstWavelenght_nm;
                const float illu_value = _illuminantSPD[idx_illu];

                assert(idx_illu < _illuminantSPD.size());

                const float interp_illu = Util::alpha(
                    wavelengths_nm[wl_idx_i],
                    wavelengths_nm[wl_idx_i + 1],
                    curr_wl_i);

                for (size_t idx_cmf = idx_cmf_start; idx_cmf < idx_cmf_end;
                     idx_cmf++) {
                    const float curr_wl_o = cmfWavelengthValue(idx_cmf);

                    const float interp_rerad = Util::alpha(
                        wavelengths_nm[wl_idx_o],
                        wavelengths_nm[wl_idx_o + 1],
                        curr_wl_o);

                    // TODO.. not sure if this is correct
                    if (wl_idx_i == wl_idx_o) {
                        normalisation_factor += illu_value * _xyzCmfs[1][idx_cmf]; // Y
                    }

                    // interpolation
                    const float bispect = (1 - interp_rerad) * ((1 - interp_illu) * value(wl_idx_i, wl_idx_o) + (interp_illu)*value(wl_idx_i + 1, wl_idx_o)) + (interp_rerad) * ((1 - interp_illu) * value(wl_idx_i, wl_idx_o + 1) + (interp_illu)*value(wl_idx_i + 1, wl_idx_o + 1));

                    const float curr_value = illu_value * bispect;

                    for (size_t c = 0; c < 3; c++) {
                        XYZ[c] += curr_value * _xyzCmfs[c][idx_cmf];
                    }
                }
            }
        }
    }

    for (size_t c = 0; c < 3; c++) {
        XYZ[c] /= normalisation_factor;
    }
}


void SpectrumConverter::spectrumToRGB(
    const std::vector<float> &wavelengths_nm,
    const float *diagonal,
    const float *reradiation,
    std::array<float, 3> &RGB) const
{
    std::array<float, 3> XYZ;
    spectrumToXYZ(wavelengths_nm, diagonal, reradiation, XYZ);

    memset(&RGB[0], 0, 3 * sizeof(float));

    // Convert to RGB using the provided matrix
    for (size_t channel = 0; channel < 3; channel++) {
        for (size_t col = 0; col < 3; col++) {
            RGB[channel] += XYZ[col] * _xyzToRgb[3 * channel + col];
        }

        // Ensure RGB values are > 0
        RGB[channel] = std::max(RGB[channel], 0.F);
    }
}


void SpectrumConverter::spectraToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *diagonal,
    const float *reradiation,
    const float *emissiveSpectrum,
    std::array<float, 3> &XYZ) const
{
    std::array<float, 3> XYZ_refl;
    std::array<float, 3> XYZ_emissive;

    spectrumToXYZ(wavelengths_nm, diagonal, reradiation, XYZ_refl);
    emissiveSpectrumToXYZ(wavelengths_nm, emissiveSpectrum, XYZ_emissive);

    for (size_t c = 0; c < 3; c++) {
        XYZ[c] = XYZ_refl[c] + XYZ_emissive[c];
    }
}


void SpectrumConverter::spectraToRGB(
    const std::vector<float> &wavelengths_nm,
    const float *diagonal,
    const float *reradiation,
    const float *emissiveSpectrum,
    std::array<float, 3> &RGB) const
{
    std::array<float, 3> XYZ;
    spectraToXYZ(
        wavelengths_nm,
        diagonal,
        reradiation,
        emissiveSpectrum,
        XYZ);

    memset(&RGB[0], 0, 3 * sizeof(float));

    // Convert to RGB using the provided matrix
    for (size_t channel = 0; channel < 3; channel++) {
        for (size_t col = 0; col < 3; col++) {
            RGB[channel] += XYZ[col] * _xyzToRgb[3 * channel + col];
        }

        // Ensure RGB values are > 0
        RGB[channel] = std::max(RGB[channel], 0.F);
    }
}


void SpectrumConverter::emissiveSpectrumToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *spectrum,
    std::array<float, 3> &XYZ) const
{
    memset(&XYZ[0], 0, 3 * sizeof(float));

    if (wavelengths_nm.size() == 0) {
        return;
    }

    const float start_wavelength = std::max(firstWavelength(), wavelengths_nm.front());
    const float end_wavelength = std::min(lastWavelength(), wavelengths_nm.back());

    // Early exit, selection out of range
    if (end_wavelength < start_wavelength) {
        return;
    }

    assert(start_wavelength <= end_wavelength);

    for (size_t idx_value = 0; idx_value < wavelengths_nm.size() - 1;
         idx_value++) {
        float wl_a = wavelengths_nm[idx_value];
        float wl_b = wavelengths_nm[idx_value + 1];

        // We have not reached yet the starting point
        if (start_wavelength > wl_b) {
            continue;
        }

        // We have finished the integration
        if (end_wavelength < wl_a) {
            break;
        }

        if (start_wavelength > wl_a) {
            wl_a = start_wavelength;
        }

        if (end_wavelength < wl_b) {
            wl_b = end_wavelength;
        }

        const size_t idx_cmf_start = cmfWavelengthIndex(wl_a);
        size_t idx_cmf_end = cmfWavelengthIndex(wl_b);

        // On last intervall we need to include the last
        // wavelength of the spectrum
        if (idx_value == wavelengths_nm.size() - 2) {
            idx_cmf_end = idx_cmf_end + 1;
        }

        for (size_t idx_cmf = idx_cmf_start; idx_cmf < idx_cmf_end;
             idx_cmf++) {
            const float curr_wl = cmfWavelengthValue(idx_cmf);
            const float curr_value = Util::interp(
                curr_wl,
                wavelengths_nm[idx_value],
                wavelengths_nm[idx_value + 1],
                spectrum[idx_value],
                spectrum[idx_value + 1]);

            for (size_t c = 0; c < 3; c++) {
                XYZ[c] += _xyzCmfs[c][idx_cmf] * curr_value;
            }
        }
    }
}


void SpectrumConverter::reflectiveSpectrumToXYZ(
    const std::vector<float> &wavelengths_nm,
    const float *spectrum,
    std::array<float, 3> &XYZ) const
{
    memset(&XYZ[0], 0, 3 * sizeof(float));

    if (wavelengths_nm.size() == 0) {
        return;
    }

    const float illuminant_last_wavelength = _illuminantFirstWavelenght_nm + _illuminantSPD.size() - 1;
    const float start_wavelength = std::max(
        std::max(_illuminantFirstWavelenght_nm, firstWavelength()),
        wavelengths_nm.front());
    const float end_wavelength = std::min(
        std::min(illuminant_last_wavelength, lastWavelength()),
        wavelengths_nm.back());

    // Early exit, selection out of range
    if (end_wavelength < start_wavelength) {
        return;
    }

    assert(start_wavelength <= end_wavelength);

    float normalisation_factor(0);

    for (size_t idx_value = 0; idx_value < wavelengths_nm.size() - 1;
         idx_value++) {
        float wl_a = wavelengths_nm[idx_value];
        float wl_b = wavelengths_nm[idx_value + 1];

        // We have not reached yet the starting point
        if (start_wavelength > wl_b) {
            continue;
        }

        // We have finished the integration
        if (end_wavelength < wl_a) {
            break;
        }

        if (start_wavelength > wl_a) {
            wl_a = start_wavelength;
        }

        if (end_wavelength < wl_b) {
            wl_b = end_wavelength;
        }

        const size_t idx_curve_start = cmfWavelengthIndex(wl_a);
        size_t idx_curve_end = cmfWavelengthIndex(wl_b);

        // On last intervall we need to include the last
        // wavelength of the spectrum
        if (idx_value == wavelengths_nm.size() - 2) {
            idx_curve_end = idx_curve_end + 1;
        }

        for (size_t idx_curve = idx_curve_start; idx_curve < idx_curve_end;
             idx_curve++) {
            const float curr_wl = cmfWavelengthValue(idx_curve);

            const size_t idx_illu_a = curr_wl - _illuminantFirstWavelenght_nm;
            assert(curr_wl >= _illuminantFirstWavelenght_nm);
            assert(idx_illu_a < _illuminantSPD.size());

            const float illu_value = _illuminantSPD[idx_illu_a];

            normalisation_factor += illu_value * _xyzCmfs[1][idx_curve]; // Y

            const float curr_value = illu_value * Util::interp(
                                                      curr_wl,
                                                      wavelengths_nm[idx_value],
                                                      wavelengths_nm[idx_value + 1],
                                                      spectrum[idx_value],
                                                      spectrum[idx_value + 1]);

            for (size_t c = 0; c < 3; c++) {
                XYZ[c] += curr_value * _xyzCmfs[c][idx_curve];
            }
        }
    }

    for (size_t c = 0; c < 3; c++) {
        XYZ[c] /= normalisation_factor;
    }
}

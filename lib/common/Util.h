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

#pragma once

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <cstdint>

/**
 * @brief Collection of utility functions
 *
 * Misc. functions. Probably would be better to separate those in
 * different groups in the future.
 */
class Util
{
public:
    template<class Float>
    static Float interp(Float x, Float x0, Float x1, Float y0, Float y1)
    {
        return lerp(y0, y1, alpha(x0, x1, x));
    }


    template<class Float>
    static Float alpha(Float x0, Float x1, Float x)
    {
        return (x - x0) / (x1 - x0);
    }


    template<class Float>
    static Float lerp(Float a, Float b, Float t)
    {
        return a + t * (b - a);
    }


    template<typename T>
    static void linspace(T min, T max, int size, std::vector<T>& array)
    {
        array.resize(size);

        T delta = (max - min) / T(size - 1);

        for (int i = 0; i < size; i++) {
            array[i] = min + i * delta;
        }
    }


    template<typename T>
    static T clamp(T value, T low, T high)
    {
        return std::max(low, std::min(high, value));
    }


    static size_t idxFromWavelengthIdx(size_t wlFrom_idx, size_t wlTo_idx)
    {
        if (wlFrom_idx < wlTo_idx) {
            return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
        } else {
            return -1;
        }
    }


    // XYZ with each component in [0..1]
    template<typename Float>
    static void xyz_to_Lab(const Float XYZ[3], Float Lab[3])
    {
        Float f[3];

        const Float epsilon = 0.008856;
        const Float kappa   = 903.3;

        // D65
        const Float coefs[3] = {0.950489, 1., 1.08840};

        for (int i = 0; i < 3; i++) {
            Float t = XYZ[i] / coefs[i];

            if (t > epsilon) {
                f[i] = std::cbrt(t);
            } else {
                f[i] = (kappa * t + 16.) / 116.;
            }
        }

        Lab[0] = 116. * f[1] - 16.;
        Lab[1] = 500. * (f[0] - f[1]);
        Lab[2] = 200. * (f[1] - f[2]);
    }


    template<typename Float>
    static Float deltaE2000(const Float Lab_1[3], const Float Lab_2[3])
    {
        const Float L1 = Lab_1[0];
        const Float a1 = Lab_1[1];
        const Float b1 = Lab_1[2];

        const Float L2 = Lab_2[0];
        const Float a2 = Lab_2[1];
        const Float b2 = Lab_2[2];

        const Float C_star_1 = std::sqrt(a1 * a1 + b1 * b1);
        const Float C_star_2 = std::sqrt(a2 * a2 + b2 * b2);
        const Float bar_C_star = (C_star_1 + C_star_2) / 2.;

        Float bar_C_start_p7 = bar_C_star * bar_C_star;   // ^2
        bar_C_start_p7 = bar_C_start_p7 * bar_C_start_p7 * bar_C_start_p7 * bar_C_star; // ^7

        const Float G = .5 * (1. - std::sqrt((bar_C_start_p7 / (bar_C_start_p7 + 6103515625.))));
        const Float a_prime_1 = (1. + G) * a1;
        const Float a_prime_2 = (1. + G) * a2;
        // std::cout << ">> " << bar_C_start_p7 / (bar_C_start_p7 + 6103515625.) << std::endl;
        const Float C_prime_1 = std::sqrt(a_prime_1 * a_prime_1 + b1 * b1);
        const Float C_prime_2 = std::sqrt(a_prime_2 * a_prime_2 + b2 * b2);

        Float h_prime_1_rad = 0;

        if (b1 != 0 || a_prime_1 != 0) {
            h_prime_1_rad = std::atan2(b1, a_prime_1);

            if (h_prime_1_rad < 0) {
                h_prime_1_rad += 2. * M_PI;
            }
        }

        Float h_prime_2_rad = 0;

        if (b2 != 0 || a_prime_2 != 0) {
            h_prime_2_rad = std::atan2(b2, a_prime_2);

            if (h_prime_2_rad < 0) {
                h_prime_2_rad += 2. * M_PI;
            }
        }

        const Float p_C_prime_12 = C_prime_1 + C_prime_2;

        Float Delta_h_prime_rad = 0;

        if (p_C_prime_12 != 0) {
            if (std::abs(h_prime_2_rad - h_prime_1_rad) <= M_PI) {
                Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad;
            } else if (h_prime_2_rad - h_prime_1_rad > M_PI) {
                Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad - 2. * M_PI;
            } else {
                Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad + 2. * M_PI;
            }
        }

        const Float Delta_L_prime = L2 - L1;
        const Float Delta_C_prime = C_prime_2 - C_prime_1;
        const Float Delta_H_prime = 2. * std::sqrt(C_prime_1 * C_prime_2) * std::sin(Delta_h_prime_rad / 2.);

        const Float bar_L_prime = (L1 + L2) / 2.;
        const Float bar_C_prime = (C_prime_1 + C_prime_2) / 2.;

        Float bar_h_prime_rad = 0;

        if (p_C_prime_12 == 0) {
            bar_h_prime_rad = h_prime_1_rad + h_prime_2_rad;
        } else if (std::abs(h_prime_1_rad - h_prime_2_rad) <= M_PI) {
            bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad) / 2.;
        } else if (h_prime_1_rad + h_prime_2_rad < 2. * M_PI) {
            bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad + 2. * M_PI) / 2.;
        } else {
            bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad - 2. * M_PI) / 2.;
        }

        const Float T = 1.
            - 0.17 * std::cos(     bar_h_prime_rad - M_PI / 6.)
            + 0.24 * std::cos(2. * bar_h_prime_rad)
            + 0.32 * std::cos(3. * bar_h_prime_rad + M_PI / 30.)
            - 0.20 * std::cos(4. * bar_h_prime_rad - 7. * M_PI / 20.);

        const Float exp_v = (bar_h_prime_rad * 180. / M_PI - 275.) / 25.;
        const Float Delta_theta_rad = M_PI / 6. * std::exp(-exp_v * exp_v);

        Float bar_C_prime_p7 = bar_C_prime * bar_C_prime;
        bar_C_prime_p7 = bar_C_prime_p7 * bar_C_prime_p7 * bar_C_prime_p7 * bar_C_prime;

        Float bar_L_prime_m50_p2 = bar_L_prime - 50.;
        bar_L_prime_m50_p2 =bar_L_prime_m50_p2 * bar_L_prime_m50_p2;

        const Float R_C = 2. * std::sqrt(bar_C_prime_p7 / (bar_C_prime_p7 + 6103515625.));
        const Float S_L = 1. + 0.015 * bar_L_prime_m50_p2 / std::sqrt(20. + bar_L_prime_m50_p2);
        const Float S_C = 1. + 0.045 * bar_C_prime;
        const Float S_H = 1. + 0.015 * bar_C_prime * T;
        const Float R_T = -std::sin(2. * Delta_theta_rad) * R_C;

        const Float K_L = 1.;
        const Float K_C = 1.;
        const Float K_H = 1.;

        const Float delta_L_r = Delta_L_prime / (K_L * S_L);
        const Float delta_C_r = Delta_C_prime / (K_C * S_C);
        const Float delta_H_r = Delta_H_prime / (K_H * S_H);

        return std::sqrt(
            delta_L_r * delta_L_r
            + delta_C_r * delta_C_r
            + delta_H_r * delta_H_r
            + R_T * delta_C_r * delta_H_r);
    }


    static void split_extension(const char* filename, std::string& base_str, std::string& extension_str);

    template<typename T, typename Q>
    static void cast_vector(const std::vector<T>& in, std::vector<Q>& out) {
        out.reserve(in.size());

        for (const T& v: in) {
            out.push_back((Q)v);
        }
    }


    template<typename T>
    static T avg_err_images(
        const std::vector<T>& reference,
        const std::vector<T>& comparison,
        size_t n_pixels,
        size_t n_bands)
    {
        T error = 0;

        for (size_t px = 0; px < n_pixels; px++) {
            T px_err = 0;

            for (size_t b = 0; b < n_bands; b++) {
                const T ref = reference [px * n_bands + b];
                const T cmp = comparison[px * n_bands + b];

                const T diff = ref - cmp;

                px_err += diff * diff;
            }

            error += std::sqrt(px_err);
        }

        return error /= (T)(n_bands * n_pixels);
    }


    template<typename T>
    static T rmse_images(
        const std::vector<T>& reference,
        const std::vector<T>& comparison,
        size_t n_pixels,
        size_t n_bands)
    {
        T error = 0;

        for (size_t px = 0; px < n_pixels; px++) {
            for (size_t b = 0; b < n_bands; b++) {
                const T ref = reference [px * n_bands + b];
                const T cmp = comparison[px * n_bands + b] ;

                const T diff = ref - cmp;

                error += diff * diff;
            }
        }

        return std::sqrt(error / (T)(n_pixels * n_bands));
    }


    template<typename T>
    static T rrmse_images(
        const std::vector<T>& reference,
        const std::vector<T>& comparison,
        size_t n_pixels,
        size_t n_bands)
    {
        T error = 0;

        for (size_t px = 0; px < n_pixels; px++) {
            T px_sum_err = 0;
            T avg = 0;

            for (size_t b = 0; b < n_bands; b++) {
                const T ref = reference [px * n_bands + b];
                const T cmp = comparison[px * n_bands + b] ;

                const T diff = ref - cmp;

                avg += ref;
                px_sum_err += diff * diff;
            }

            // TODO: Double check
            px_sum_err /= (T)n_bands;
            avg /= (T)n_bands;

            if (avg > 0) {
                error += std::sqrt(px_sum_err) / avg;
            }
        }

        return error / (T)n_pixels;
    }


    template<typename T>
    static T max_error_images(
        const std::vector<T>& reference,
        const std::vector<T>& comparison,
        size_t n_pixels,
        size_t n_bands)
    {
        T max_error = 0;

        for (size_t px = 0; px < n_pixels; px++) {
            for (size_t b = 0; b < n_bands; b++) {
                const T ref = reference [px * n_bands + b];
                const T cmp = comparison[px * n_bands + b] ;

                const T diff = ref - cmp;

                max_error = std::max(max_error, std::abs(diff));
            }
        }

        return max_error;
    }
    // // TODO: Check this
    // template<typename T>
    // static T error_images(
    //     const std::vector<T>& reference,
    //     const std::vector<T>& comparison,
    //     size_t n_pixels,
    //     size_t n_bands)
    // {
    //     T error = 0;

    //     for (size_t px = 0; px < n_pixels; px++) {
    //         T px_sum_err = 0;
    //         // T avg = 0;

    //         for (size_t b = 0; b < n_bands; b++) {
    //             const T ref = reference [px * n_bands + b];
    //             const T cmp = comparison[px * n_bands + b] ;

    //             const T diff = ref - cmp;

    //             // avg += ref;
    //             px_sum_err += diff * diff;
    //         }

    //         // avg /= (T)n_bands;

    //         // if (avg > 0) {
    //         //     error += std::sqrt(px_sum_err) / avg;
    //         // }

    //         error += std::sqrt(px_sum_err);
    //     }

    //     return error / (T)n_pixels;
    // }


    template<typename T>
    static T quantize_dequantize(T src, size_t n_bits)
    {
        const uint64_t v = (1 << n_bits) - 1;

        return std::round(src * (T)v) / (T)v;
    }


    static bool ends_with(std::string const & value, std::string const & ending)
    {
        if (ending.size() > value.size()) return false;
        return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    }
};

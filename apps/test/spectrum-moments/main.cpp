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

#include <iostream>
#include <fstream>
#include <vector>
#include <cstddef>
#include <cassert>

#include "moments.h"

template<typename T>
void linspace(T min, T max, size_t size, std::vector<T>& array)
{
    array.resize(size);

    T delta = (max - min) / T(size - 1);

    for (size_t i = 0; i < size; i++) {
        array[i] = min + i * delta;
    }
}


void print_array(const std::vector<float>& array)
{
    std::cout << "[";

    for (size_t i = 0; i < array.size(); i++) {
        printf(" %.8f ", array[i]);
    }

    std::cout << "]" << std::endl;
}


void print_array(const std::vector<double>& array)
{
    std::cout << "[";

    for (size_t i = 0; i < array.size(); i++) {
        printf(" %.8lf ", array[i]);
    }

    std::cout << "]" << std::endl;
}


double get_error(const std::vector<double>& arr1, const std::vector<double>& arr2)
{
    assert(arr1.size() == arr2.size());

    double err = 0;

    for (size_t i = 0; i < arr1.size(); i++) {
        const double curr_err = arr1[i] - arr2[i];
        err += curr_err * curr_err;
    }

    return err;
}


void dump_array(const char* filename, const std::vector<double>& array)
{
    FILE* f = std::fopen(filename, "wb");

    if (!f) {
        std::cerr << "Could not create file for dumping array" << std::endl;
        return;
    }

    const uint32_t array_sz = array.size();
    std::fwrite(&array_sz, sizeof(uint32_t), 1, f);
    std::fwrite(array.data(), sizeof(double), array_sz, f);

    std::fclose(f);
}


void gnuplot_export(const char* filename, const std::vector<double>& array)
{
    std::ofstream f_out(filename, std::ios::out);

    if (!f_out.is_open()) {
        std::cerr << "Could not open file for writing." << std::endl;
        return;
    }

    for (size_t i = 0; i < array.size(); i++) {
        f_out << i << " " << array[i] << std::endl;
    }

    f_out.close();
}


int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    // ------------------------------------------------------------------------
    // Data
    // ------------------------------------------------------------------------

    std::vector<double> reflectance_wavelength;
    linspace(400., 700., 61, reflectance_wavelength);

    std::vector<double> reflectance = {
        0.060000000, 0.061000000, 0.061000000, 0.061000000, 0.062000000, 0.063000000, 0.064000000, 0.066000000,
        0.068000000, 0.071000000, 0.075000000, 0.079000000, 0.085000000, 0.093000000, 0.104000000, 0.118000000,
        0.135000000, 0.157000000, 0.185000000, 0.221000000, 0.269000000, 0.326000000, 0.384000000, 0.440000000,
        0.484000000, 0.516000000, 0.534000000, 0.542000000, 0.545000000, 0.541000000, 0.533000000, 0.524000000,
        0.513000000, 0.501000000, 0.487000000, 0.472000000, 0.454000000, 0.436000000, 0.416000000, 0.394000000,
        0.374000000, 0.358000000, 0.346000000, 0.337000000, 0.331000000, 0.328000000, 0.325000000, 0.322000000,
        0.320000000, 0.319000000, 0.319000000, 0.320000000, 0.324000000, 0.330000000, 0.337000000, 0.345000000,
        0.354000000, 0.362000000, 0.368000000, 0.375000000, 0.379000000
    };

    std::vector<double> emission_wavelength;
    linspace(312.5, 811.5, 999, emission_wavelength);

    std::vector<double> emission = {
        0.03183, 0.02949, 0.02075, 0.00807, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
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
        0.01058, 0.02523, 0.04864
    };

    const double FL3_13_start = 380.0;
    const double FL3_13_end = 780.0;
    const int FL3_13_n_samples = 81;
    const std::vector<double> FL3_13 = {2.23, 2.92, 3.91, 4.55, 5.46, 77.4, 11.25, 7.69, 8.29, 8.98, 10.01, 204.45, 13.75, 16.88, 21.73, 27.96, 34.92, 41.96, 48.62, 54.33, 59.49, 67.91, 70.01, 66.4, 62.07, 56.95, 52.7, 48.54, 44.8, 41.75, 39.77, 40.5, 59.27, 184.09, 59.06, 49.95, 50.9, 54.51, 58.33, 77.49, 85.78, 76.2, 78.73, 78.95, 81.48, 84.57, 87.75, 89.56, 91.36, 89.0, 83.67, 78.26, 73.19, 67.61, 61.42, 55.49, 49.78, 44.46, 39.13, 34.45, 30.28, 26.37, 23.88, 20.1, 17.4, 15.29, 13.62, 11.68, 10.31, 9.11, 8.03, 7.13, 6.31, 5.67, 5.11, 4.55, 9.06, 3.74, 4.04, 3.14, 2.75};

    std::vector<double> FL3_13_wavelengths;
    linspace(FL3_13_start, FL3_13_end, FL3_13_n_samples, FL3_13_wavelengths);

    std::vector<double> phases;
    std::vector<double> moments;
    std::vector<double> signal_back;

    wavelengths_to_phases(FL3_13_wavelengths, phases);
    compute_moments(phases, FL3_13, phases.size(), moments);
    compute_density(phases, moments, signal_back);

    gnuplot_export("signal.dat", FL3_13);
    gnuplot_export("reconst.dat", signal_back);
    // dump_array("signal.dat", FL3_13);
    // dump_array("reconst.dat", signal_back);


    // // ------------------------------------------------------------------------
    // // Test with reflectance
    // // ------------------------------------------------------------------------

    // const int n_moments_reflectance = 4;

    // std::vector<float> reflectance_phases;

    // std::vector<float> reflectance_moments;
    // std::vector<float> reflectance_compressed;

    // std::vector<float> reflectance_inflated;
    // std::vector<float> reflectance_density;

    // wavelengths_to_phases(reflectance_wavelength, reflectance_phases);

    // compute_moments(
    //     reflectance_phases,
    //     reflectance,
    //     n_moments_reflectance,
    //     reflectance_moments
    // );

    // bounded_compress_moments(
    //     reflectance_moments,
    //     reflectance_compressed
    // );

    // bounded_decompress_moments(
    //     reflectance_compressed,
    //     reflectance_inflated
    // );

    // bounded_compute_density_lagrange(
    //     reflectance_phases,
    //     reflectance_inflated,
    //     reflectance_density
    // );

    // std::ofstream out_reflectance("reflectance.dat");

    // for (size_t i = 0; i < reflectance_wavelength.size(); i++) {
    //     out_reflectance << reflectance_wavelength[i] << " "
    //                     << reflectance[i] << " "
    //                     << reflectance_density[i] << std::endl;
    // }

    // ------------------------------------------------------------------------
    // Test with emission
    // ------------------------------------------------------------------------

    // const int n_moments_emission = FL3_13_n_samples;

    // double max_signal = 0;
    // for (const double& v: emission) {
    //     max_signal = std::max(max_signal, v);
    // }

    // for (double& v: emission) {
    //     v /= max_signal;
    // }

    // std::vector<double> emission_phases;

    // std::vector<double> emission_moments;
    // std::vector<double> emission_compressed_bounded;
    // std::vector<double> emission_compressed_unbounded;

    // std::vector<double> emission_inflated_unbounded;
    // std::vector<double> emission_inflated_bounded;

    // std::vector<double> emission_density;

    // wavelengths_to_phases(emission_wavelength, emission_phases);

    // compute_moments(
    //     emission_phases,
    //     emission,
    //     n_moments_emission,
    //     emission_moments
    // );

    // // Unbounded
    // unbounded_compress_moments(
    //     emission_moments,
    //     emission_compressed_unbounded
    // );

    // unbounded_decompress_moments(
    //     emission_compressed_unbounded,
    //     emission_inflated_unbounded
    // );

    // // Bounded
    // bounded_compress_moments(
    //     emission_moments,
    //     emission_compressed_bounded
    // );

    // bounded_decompress_moments(
    //     emission_compressed_bounded,
    //     emission_inflated_bounded
    // );

    // compute_density(
    //     emission_phases,
    //     emission_inflated_bounded,
    //     emission_density
    // );

    // std::cout << "-- Moments --" << std::endl;
    // print_array(emission_moments);
    // std::cout << "-- Moments compressed unbounded --" << std::endl;
    // print_array(emission_compressed_unbounded);
    // std::cout << "-- Moments compressed bounded--" << std::endl;
    // print_array(emission_compressed_unbounded);

    // std::cout << "------------------------------------" << std::endl;
    // std::cout << "--   Comparison forward reverse   --" << std::endl;
    // std::cout << "------------------------------------" << std::endl;
    // std::cout << "  bounded original: " << get_error(emission_moments, emission_inflated_bounded) << std::endl;
    // std::cout << "unbounded original: " << get_error(emission_moments, emission_inflated_unbounded) << std::endl;

    // dump_array("phases.dat", emission_phases);
    // dump_array("moments.dat", emission_moments);
    // dump_array("c_moments_b.dat", emission_compressed_bounded);
    // dump_array("c_moments_ub.dat", emission_compressed_unbounded);

    // print_array(emission_moments);

    // std::ofstream out_emission("emission.dat");

    // for (size_t i = 0; i < emission_wavelength.size(); i++) {
    //     out_emission << emission_wavelength[i] << " "
    //                  << emission[i] << " "
    //                  << emission_density[i] << std::endl;
    // }

    return 0;
}

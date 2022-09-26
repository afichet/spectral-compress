#include <iostream>

#include "data.h"
#include <moments.h>

#include <cassert>
#include <vector>
#include <Eigen/Core>

void print_array(const std::vector<float>& array)
{
    std::cout << "[";

    for (size_t i = 0; i < array.size(); i++) {
        printf(" %.8f ", array[i]);
    }

    std::cout << "]" << std::endl;
}

void print_array_diff(
    const std::vector<float>& array_1,
    const std::vector<float>& array_2)
{
    assert(array_1.size() == array_2.size());

    std::cout << "[";

    for (size_t i = 0; i < array_1.size(); i++) {
        printf(" %.8f ", array_1[i] - array_2[i]);
    }

    std::cout << "]" << std::endl;
}


int main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;

    // Reflectance spectrum
    std::vector<float> reflectance_wavelength;
    linspace(400.f, 700.f, 61, reflectance_wavelength);

    std::vector<float> reflectance_phases;
    std::vector<float> reflectance_moments;
    std::vector<float> reflectance_compressed;

    std::cout << "-- PHASES --" << std::endl;
    wavelengths_to_phases(reflectance_wavelength, reflectance_phases);
    print_array(reflectance_phases);

    std::cout << "-- MOMENTS --" << std::endl;
    compute_moments(reflectance_phases, reflectance, reflectance_phases.size(), reflectance_moments);
    print_array(reflectance_moments);

    std::cout << "-- COMPRESSED BOUNDED --" << std::endl;
    compress_bounded_moments(reflectance_moments, reflectance_compressed);
    print_array(reflectance_compressed);


    std::vector<float> b_reflectance_moments;
    std::vector<float> b_reflectance(reflectance_moments.size());
    std::vector<float> basis(reflectance_moments.size() * reflectance_moments.size());

    std::cout << "-- BACKWARD MOMENTS --" << std::endl;
    decompress_bounded_moments(reflectance_compressed, b_reflectance_moments);
    print_array_diff(reflectance_moments, b_reflectance_moments);

    std::cout << "-- BACKWARD SIGNAL --" << std::endl;
    compute_basis_moments_to_signal(reflectance_phases, basis);
    Eigen::Map<Eigen::MatrixXf> t_mat(basis.data(), reflectance_moments.size(), reflectance_moments.size());
    Eigen::Map<Eigen::VectorXf> signal(b_reflectance.data(), reflectance_moments.size());

    const Eigen::Map<Eigen::VectorXf> moments(b_reflectance_moments.data(), reflectance_moments.size());
    signal = t_mat * moments;
    print_array_diff(reflectance, b_reflectance);


    // Emission spectrum
    std::vector<float> emission_wavelength;
    linspace(312.5f, 811.5f, 999, emission_wavelength);


    return 0;
}

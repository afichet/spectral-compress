#include <iostream>

#include "data.h"
#include <moments.h>

#include <cassert>
#include <vector>
#include <Eigen/Core>

void print_array(const std::vector<double>& array)
{
    std::cout << "[";

    for (size_t i = 0; i < array.size(); i++) {
        printf(" %.8f ", array[i]);
    }

    std::cout << "]" << std::endl;
}

void print_array_diff(
    const std::vector<double>& array_1,
    const std::vector<double>& array_2)
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
    std::vector<double> reflectance_wavelength;
    linspace(400.f, 700.f, 61, reflectance_wavelength);

    std::vector<double> reflectance_phases;
    std::vector<double> reflectance_moments;
    std::vector<double> reflectance_compressed;

    std::cout << "-- PHASES --" << std::endl;
    wavelengths_to_phases(reflectance_wavelength, reflectance_phases);
    print_array(reflectance_phases);

    std::cout << "-- MOMENTS --" << std::endl;
    compute_moments(reflectance_phases, reflectance, reflectance_phases.size(), reflectance_moments);
    print_array(reflectance_moments);

    std::cout << "-- COMPRESSED BOUNDED --" << std::endl;
    unbounded_to_bounded_compress_moments(reflectance_moments, reflectance_compressed);
    print_array(reflectance_compressed);


    std::vector<double> b_reflectance_moments;
    std::vector<double> b_reflectance(reflectance_moments.size());
    std::vector<double> basis(reflectance_moments.size() * reflectance_moments.size());

    std::cout << "-- BACKWARD MOMENTS --" << std::endl;
    unbounded_to_bounded_decompress_moments(reflectance_compressed, b_reflectance_moments);
    print_array_diff(reflectance_moments, b_reflectance_moments);

    std::cout << "-- BACKWARD SIGNAL --" << std::endl;
    compute_basis_moments_to_signal(reflectance_phases, basis);
    
    Eigen::Map<Eigen::MatrixXd> t_mat(basis.data(), reflectance_moments.size(), reflectance_moments.size());
    Eigen::Map<Eigen::VectorXd> signal(b_reflectance.data(), reflectance_moments.size());
    const Eigen::Map<Eigen::VectorXd> moments(b_reflectance_moments.data(), reflectance_moments.size());

    signal = t_mat * moments;
    print_array_diff(reflectance, b_reflectance);

    std::cout << " ================================================= " << std::endl;
    std::cout << " ================================================= " << std::endl;
    std::cout << " ================================================= " << std::endl;

    // Emission spectrum
    std::vector<double> emission_wavelength;
    linspace(312.5f, 811.5f, 999, emission_wavelength);

    std::vector<double> emission_phases;
    std::vector<double> emission_moments;
    std::vector<double> emission_compressed;

    std::cout << "-- PHASES --" << std::endl;
    wavelengths_to_phases(emission_wavelength, emission_phases);
    print_array(emission_phases);

    std::cout << "-- MOMENTS --" << std::endl;
    compute_moments(emission_phases, emission, emission_phases.size(), emission_moments);
    print_array(emission_moments);

    std::cout << "-- COMPRESSED U->B --" << std::endl;
    // unbounded_compress_moments(emission_moments, emission_compressed);
    unbounded_to_bounded_compress_moments(emission_moments, emission_compressed);
    print_array(emission_compressed);

    std::vector<double> b_emission_moments;
    std::vector<double> b_emission(emission_phases.size());
    std::vector<double> basis_emission(emission_phases.size() * emission_phases.size());

    std::cout << "-- BACKWARD MOMENTS --" << std::endl;
    // unbounded_decompress_moments(emission_compressed, b_emission_moments);
    unbounded_to_bounded_decompress_moments(emission_compressed, b_emission_moments);
    print_array_diff(emission_moments, b_emission_moments);

    std::cout << "-- BACKWARD SIGNAL --" << std::endl;
    compute_basis_moments_to_signal(emission_phases, basis_emission);

    Eigen::Map<Eigen::MatrixXd> e_t_mat(basis_emission.data(), emission_moments.size(), emission_moments.size());
    Eigen::Map<Eigen::VectorXd> e_signal(b_emission.data(), emission_moments.size());
    const Eigen::Map<Eigen::VectorXd> e_moments(b_emission_moments.data(), emission_moments.size());

    e_signal = e_t_mat * e_moments;
    print_array_diff(emission, b_emission);


    return 0;
}

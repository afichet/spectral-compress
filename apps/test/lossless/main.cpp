#include <iostream>
#include <vector>

#include <Eigen/Core>

#include <moments.h>


void linspace(float min, float max, size_t size, std::vector<float>& array)
{
    array.resize(size);

    float delta = (max - min) / float(size - 1);

    for (size_t i = 0; i < size; i++) {
        array[i] = min + i * delta;
    }
}


int main(int argc, char* argv[])
{
    std::vector<float> reflectance = {
        0.060, 0.062, 0.068, 0.085,
        0.135, 0.269, 0.484, 0.545, 
        0.513, 0.454, 0.374, 0.331, 
        0.320, 0.324, 0.354, 0.379
    };

    std::vector<float> reflectance_wavelength;
    linspace(400.f, 700.f, reflectance.size(), reflectance_wavelength);

    std::vector<float> reflectance_phases;
    wavelengths_to_phases(reflectance_wavelength, reflectance_phases);

    const size_t n_moments = reflectance.size();

    // Directly compute the transformation
    std::vector<float> moments;
    compute_moments(reflectance_phases, reflectance, n_moments - 1, moments);
    
    Eigen::Map<Eigen::RowVectorXf> builtin_moments(moments.data(), moments.size());

    std::cout << "Built in" << std::endl;
    std::cout << builtin_moments << std::endl;

    // Construct the transformation matrix    
    std::vector<float> std_signal_moments(n_moments * n_moments);
    compute_basis_signal_to_moments(reflectance_phases, std_signal_moments);

    Eigen::Map<Eigen::MatrixXf> transform_mat(std_signal_moments.data(), n_moments, n_moments);    
    Eigen::Map<Eigen::VectorXf> refl(reflectance.data(), reflectance.size());
    
    Eigen::VectorXf eigen_moments = transform_mat * refl;

    std::cout << "Eigen: " << std::endl;
    std::cout << eigen_moments.transpose() << std::endl;

    std::vector<float> std_moments_signal(n_moments * n_moments);
    compute_basis_moments_to_signal(reflectance_phases, std_moments_signal);

    Eigen::Map<Eigen::MatrixXf> transform_inv(std_moments_signal.data(), n_moments, n_moments);

    Eigen::VectorXf eigen_signal = transform_inv * eigen_moments;
    
    std::cout << "Signal retrieved:" << std::endl;
    std::cout << eigen_signal.transpose() << std::endl;

    return 0;
}
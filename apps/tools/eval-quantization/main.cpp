#include <iostream>

#include <EXRSpectralImage.h>
#include <quantization.h>
#include <sstream>
#include <fstream>

int main(int argc, char* argv[]) 
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral exr> <output_prefix>" << std::endl;
        return 0;
    }

    // const int n_bits_1 = 10;
    const char* filename = argv[1];
    const std::string output_prefix = argv[2];

    EXRSpectralImage image(filename);

    const int n_pixels = image.width() * image.height();

    const int n_bits[] = {12, 10, 8};

    for (SpectralFramebuffer* fb: image.getSpectralFramebuffers()) {
        for (int bits: n_bits) {
            const int n_bands = fb->wavelengths_nm.size();
            std::vector<int> quantization_curve_b;

            std::vector<double> wavelengths(n_bands);
            std::vector<double> image_data(n_pixels * n_bands);
            
            // Cast to double
            for (int i = 0; i < n_bands; i++) {
                wavelengths[i] = (double)fb->wavelengths_nm[i];
            }

            for (int i = 0; i < n_pixels * n_bands; i++) {
                image_data[i] = (double)fb->image_data[i];
            }

            int start_n_bits = fb->pixel_type == PixelType::HALF ? 16 : 32;

            double err_utb = bounded_compute_quantization_curve(
                wavelengths,
                image_data,
                n_pixels, n_bands, 
                bits, 
                quantization_curve_b,
                start_n_bits
            );

            // Save data
            std::stringstream output_file;
            output_file << output_prefix << "_b_" << bits << ".txt";

            std::ofstream out_utb(output_file.str());

            for (int i = 0; i < quantization_curve_b.size(); i++) {
                out_utb << quantization_curve_b[i] << " ";
            }

            out_utb << std::endl << err_utb;
        }
    }


    return 0;
}
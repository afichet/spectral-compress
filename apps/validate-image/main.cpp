#include <iostream>
#include <EXRSpectralImage.h>

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_image>" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];

    const SEXR::EXRSpectralImage image_in(filename_in);

    const size_t width     = image_in.width();
    const size_t height    = image_in.height();
    const size_t n_bands   = image_in.nSpectralBands();

    size_t n_nans     = 0;
    size_t n_inf      = 0;
    size_t n_neg      = 0;
    size_t n_all_zero = 0;

    const float* framebuffer = &image_in.emissive(0, 0, 0, 0);

    for (size_t i = 0; i < width * height; i++) {
        bool all_zeros = true;

        for (size_t b = 0; b < n_bands; b++) {
            const float v = framebuffer[n_bands * i + b];

            if (std::isinf(v)) {
                ++n_inf;
            } else if (std::isnan(v)) {
                ++n_nans;
            } else if (v < 0) {
                ++n_neg;
            }
            else if (v > 0) {
                all_zeros = false;
            }
        }

        if (all_zeros) {
            ++n_all_zero;
        }
    }

    std::cout << "Summary:" << std::endl;
    std::cout << "        NaNs: " << n_nans << std::endl;
    std::cout << "         Inf: " << n_inf << std::endl;
    std::cout << "         Neg: " << n_neg << std::endl;
    std::cout << "Zero spectra: " << n_all_zero << std::endl;



    return 0;
}
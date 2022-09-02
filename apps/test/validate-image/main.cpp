#include <iostream>
#include <cmath>
#include <EXRImage.h>

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_image>" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];

    const EXRImage exr_in(filename_in);

    const size_t width     = exr_in.width();
    const size_t height    = exr_in.height();

    size_t n_nans     = 0;
    size_t n_inf      = 0;
    size_t n_neg      = 0;

    const std::vector<EXRFramebuffer*>& framebuffers = exr_in.getFramebuffersConst();

    for (const EXRFramebuffer* fb: framebuffers) {
        for (size_t i = 0; i < width * height; i++) {
            const float v = fb->getPixelDataConst()[i];

            if (std::isinf(v)) {
                ++n_inf;
            } else if (std::isnan(v)) {
                ++n_nans;
            } else if (v < 0) {
                ++n_neg;
            }
        }
    }

    std::cout << "Summary:" << std::endl;
    std::cout << "        NaNs: " << n_nans << std::endl;
    std::cout << "         Inf: " << n_inf << std::endl;
    std::cout << "         Neg: " << n_neg << std::endl;


    return 0;
}
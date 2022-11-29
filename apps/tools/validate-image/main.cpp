#include <iostream>
#include <cmath>
#include <EXRImage.h>
#include <lodepng.h>
#include <cstdint>
#include <vector>
#include <string>

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_image> [<mask_err_png>]" << std::endl;

        exit(0);
    }

    const char* filename_in  = argv[1];

    bool write_mask = false;
    char* filename_mask = NULL;

    if (argc >= 3) {
        write_mask = true;
        filename_mask = argv[2];
    }

    const EXRImage exr_in(filename_in);

    const size_t width     = exr_in.width();
    const size_t height    = exr_in.height();

    std::vector<uint8_t> mask_err(width * height * 4);
    memset(mask_err.data(), 0, width * height * 4);

    size_t n_nans     = 0;
    size_t n_inf      = 0;
    size_t n_neg      = 0;

    const std::vector<EXRFramebuffer*>& framebuffers = exr_in.getFramebuffersConst();

    for (const EXRFramebuffer* fb: framebuffers) {
        // Ignore RGB layers
        const std::string name = fb->getName();
        if (name.find("R") && name.find("G") && name.find("B")) {
            for (size_t i = 0; i < width * height; i++) {
                const float v = fb->getPixelDataConst()[i];

                if (std::isinf(v)) {
                    ++n_inf;
                    mask_err[4 * i + 0] = 255;
                } else if (std::isnan(v)) {
                    ++n_nans;
                    mask_err[4 * i + 1] = 255;
                } else if (v < 0) {
                    ++n_neg;
                    mask_err[4 * i + 2] = 255;
                }

                mask_err[4 * i + 3] = 255;
            }
        }
    }

    std::cout << "Summary:" << std::endl;
    std::cout << "        NaNs: " << n_nans << std::endl;
    std::cout << "         Inf: " << n_inf << std::endl;
    std::cout << "         Neg: " << n_neg << std::endl;

    if (write_mask) {
        lodepng::encode(filename_mask, mask_err.data(), width, height);
    }

    return 0;
}

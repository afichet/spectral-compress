#include <EXRSpectralImage.h>
#include "macbeth_data.h"

#include <iostream>
#include <vector>
#include <cstring>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <EXR image out>" << std::endl;

        exit(0);
    }

    size_t width = 6;
    size_t height = 4;
    size_t n_bands = 36;

    EXRSpectralImage image_out(width, height);

    std::vector<float> wavelengths(n_bands);
    std::memcpy(wavelengths.data(), macbeth_wavelengths, n_bands * sizeof(float));

    std::vector<float> framebuffer(width * height * n_bands);

    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            const size_t patch_idx =  y * width + x;

            std::memcpy(
                &framebuffer[n_bands * patch_idx],
                macbeth_patches[patch_idx],
                n_bands * sizeof(float)
            );
        }
    }

    image_out.appendSpectralFramebuffer(wavelengths, framebuffer, "T");

    image_out.write(argv[1]);

    return 0;
}
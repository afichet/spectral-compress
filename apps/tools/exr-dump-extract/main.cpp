#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <lodepng.h>
#include <EXRSpectralImage.h>

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <exr-dump>" << std::endl;

        return 0;
    }

    const float exposure = 5.f;
    const float gamma = 1.f/2.2f;

    EXRSpectralImage* img = EXRSpectralImage::read_dump(argv[1]);

    for(const SpectralFramebuffer* fb: img->getSpectralFramebuffers()) {
        for (size_t i = 0; i < fb->wavelengths_nm.size(); i++) {
            std::stringstream filename;
            filename << fb->root_name << "." << fb->wavelengths_nm[i] << ".png";

            std::vector<uint8_t> out_fb(4 * img->width() * img->height());

            for (size_t p = 0; p < img->width() * img->height(); p++) {
                const float pixel = fb->image_data[p * fb->wavelengths_nm.size() + i];
                const float value = std::max(0.f, std::min(1.f, std::pow(std::exp2(exposure) * pixel, gamma)));
                
                out_fb[4 * p + 0] = 255 * value;
                out_fb[4 * p + 1] = 255 * value;
                out_fb[4 * p + 2] = 255 * value;
                out_fb[4 * p + 3] = 255;
            }

            lodepng::encode(filename.str(), out_fb.data(), img->width(), img->height());
        }
    }

    delete img;

    return 0;
}
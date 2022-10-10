#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include <lodepng.h>
#include <JXLImage.h>

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <jxl-dump>" << std::endl;

        return 0;
    }

    const float exposure = 0;//5.f;
    const float gamma = 1.f/2.2f;

    JXLImage* img = JXLImage::read_dump(argv[1]);
    int i = 0;

    for(const JXLFramebuffer* fb: img->getFramebuffers()) {
        std::stringstream filename;

        if (fb->getName()) {
            ++i;
            filename << i; // fb->getName();
        } else {
            filename << "no_name";
        }

        filename << ".png";

        std::vector<uint8_t> out_fb(4 * img->width() * img->height());
        const std::vector<float>& in_fb = fb->getPixelDataConst();

        for (size_t p = 0; p < img->width() * img->height(); p++) {
            const float pixel = in_fb[p];
            const float value = std::max(0.f, std::min(1.f, std::pow(std::exp2(exposure) * pixel, gamma)));
            
            out_fb[4 * p + 0] = 255 * value;
            out_fb[4 * p + 1] = 255 * value;
            out_fb[4 * p + 2] = 255 * value;
            out_fb[4 * p + 3] = 255;
        }

        lodepng::encode(filename.str(), out_fb.data(), img->width(), img->height());
    }

    delete img;

    return 0;
}
#include <iostream>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <vector>

#include <lodepng.h>

#include <JXLImage.h>

void export_float_gray_to_png(
    const std::vector<float>& framebuffer,
    uint32_t width, uint32_t height,
    const char* filename)
{
    std::vector<uint8_t> export_fb(framebuffer.size() * 4);

    for (size_t i = 0; i < framebuffer.size(); i++) {
        const uint8_t value = 255.f * std::max(0.f, std::min(1.f, framebuffer[i]));

        for (size_t c = 0; c < 3; c++) {
            export_fb[4 * i + c] = value;
        }

        export_fb[4 * i + 3] = 255;
    }

    lodepng::encode(filename, export_fb.data(), width, height);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <jxl_in>" << std::endl;
        exit(0);
    }

    JXLImage image_in(argv[1]);

    const uint32_t width  = image_in.width();
    const uint32_t height = image_in.height();

    // Export the main framebuffer
    std::vector<float> main_fb(width * height);
    std::memcpy(main_fb.data(), image_in.getFramebufferData(0).data(), width * height * sizeof(float));

    float min_v = main_fb[0], max_v = main_fb[0];

    // Rescaling
    for (uint32_t i = 0; i < main_fb.size(); i++) {
        min_v = std::min(min_v, main_fb[i]);
        max_v = std::max(max_v, main_fb[i]);
    }

    for (uint32_t i = 0; i < main_fb.size(); i++) {
        main_fb[i] = 1000 * main_fb[i] * (max_v - min_v) + min_v;
    }

    std::stringstream ss;
    ss << argv[1] << "_0.png";

    export_float_gray_to_png(main_fb, width, height, ss.str().c_str());


    for (uint32_t layer = 1; layer < image_in.n_framebuffers(); layer++) {
        std::stringstream ss;
        ss << argv[1] << "_" << layer << ".png";

        export_float_gray_to_png(
            image_in.getFramebufferData(layer),
            width, height,
            ss.str().c_str());
    }

    return 0;
}

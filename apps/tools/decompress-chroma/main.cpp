#include <iostream>

#include <JXLImage.h>
#include <EXRSpectralImage.h>
#include <Util.h>
#include <moments_image.h>
#include <moments.h>


void decompress_spectral_framebuffer(
    const SGEGSpectralGroup& sg,
    const std::vector<std::vector<float>>& compressed_moments,
    std::vector<float>& spectral_framebuffer)
{
    const size_t n_pixels = compressed_moments[0].size();

    size_t n_moments = 0;

    switch (sg.method) {
        case LINEAR:
        case LINAVG:
        case BOUNDED:
        case UNBOUNDED:
        case UNBOUNDED_TO_BOUNDED:
            assert(sg.layer_indices.size() == (sg.mins.size() + 1));
            assert(sg.layer_indices.size() == (sg.maxs.size() + 1));
            n_moments = compressed_moments.size();
            break;
        case UPPERBOUND:
        case TWOBOUNDS:
            assert(sg.layer_indices.size() - 1 == (sg.mins.size() + 1));
            assert(sg.layer_indices.size() - 1 == (sg.maxs.size() + 1));
            n_moments = compressed_moments.size() - 1;
            break;
    }

    std::vector<double> wavelengths_d;
    std::vector<double> mins_d, maxs_d;
    std::vector<double> spectral_framebuffer_d;
    std::vector<double> compressed_moments_d(n_moments * n_pixels);
    std::vector<uint8_t> relative_scales;

    Util::cast_vector(sg.wavelengths, wavelengths_d);
    Util::cast_vector(sg.mins, mins_d);
    Util::cast_vector(sg.maxs, maxs_d);

    for (size_t px = 0; px < n_pixels; px++) {
        for (size_t m = 0; m < n_moments; m++) {
            compressed_moments_d[px * n_moments + m] = compressed_moments[m][px];
        }

        compressed_moments_d[px * n_moments] = 0.5;
    }

    if (sg.method == UPPERBOUND || sg.method == TWOBOUNDS) {
        relative_scales.resize(n_pixels);

        for (size_t px = 0; px < n_pixels; px++) {
            relative_scales[px] = std::numeric_limits<uint8_t>::max() * compressed_moments[n_moments][px];
        }
    }

    double timing;

    decompress_spectral_image(
        LINEAR,
        wavelengths_d,
        compressed_moments_d,
        mins_d, maxs_d,
        relative_scales,
        sg.global_min, sg.global_max,
        n_pixels,
        n_moments,
        spectral_framebuffer_d,
        timing
    );

    Util::cast_vector(spectral_framebuffer_d, spectral_framebuffer);
}


Imf::PixelType quantization_to_exr(size_t n_bits, size_t n_exponent_bits) {
    if (n_exponent_bits > 0) {
        if (n_bits > 16) {
            return Imf::PixelType::FLOAT;
        } else {
            return Imf::PixelType::HALF;
        }
    } else {
        return Imf::PixelType::UINT;
    }
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <jxl_in> <exr_out>" << std::endl;
        exit(0);
    }

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    const JXLImage jxl_image(filename_in);

    const SGEGBox box   = jxl_image.getBox();
    const size_t width  = jxl_image.width();
    const size_t height = jxl_image.height();


    EXRSpectralImage exr_out(width, height);

    exr_out.setAttributesData(box.exr_attributes);

    for (const SGEGSpectralGroup& sg: box.spectral_groups) {
        std::string root_name = sg.root_name.data();
        const size_t n_moments = sg.layer_indices.size();

        assert(sg.mins.size() == sg.maxs.size());

        std::vector<std::vector<float>> moments(n_moments);
        std::vector<float> spectral_framebuffer;

        for (size_t m = 0; m < n_moments; m++) {
            moments[m] = jxl_image.getFramebufferDataConst(sg.layer_indices[m]);
        }

        decompress_spectral_framebuffer(
            sg,
            moments,
            spectral_framebuffer
        );

        const JXLFramebuffer* main_fb   = jxl_image.getFramebuffer(sg.layer_indices[0]);
        const size_t n_bits             = main_fb->getBitsPerSample();
        const size_t n_exponent_bits    = main_fb->getExponentBitsPerSample();
        const Imf::PixelType pixel_type = quantization_to_exr(n_bits, n_exponent_bits);

        exr_out.appendSpectralFramebuffer(
            sg.wavelengths,
            spectral_framebuffer,
            root_name,
            pixel_type
        );
    }

    for (const SGEGGrayGroup& gg: box.gray_groups) {
        const JXLFramebuffer* main_fb   = jxl_image.getFramebuffer(gg.layer_index);
        const size_t n_bits             = main_fb->getBitsPerSample();
        const size_t n_exponent_bits    = main_fb->getExponentBitsPerSample();
        const Imf::PixelType pixel_type = quantization_to_exr(n_bits, n_exponent_bits);

        exr_out.appendExtraFramebuffer(
            jxl_image.getFramebufferDataConst(gg.layer_index),
            gg.layer_name.data(),
            pixel_type
        );
    }

    exr_out.write(filename_out);


    return 0;
}

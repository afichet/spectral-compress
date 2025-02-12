/**
 * Copyright 2022 - 2023 Alban Fichet
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 *   3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <limits>

#include <cstdint>
#include <cmath>
#include <cassert>

#include <tclap/CmdLine.h>

#include <png.h>

#include <SpectrumConverter.h>

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfFrameBuffer.h>
#include <ImfHeader.h>

// TODO: accept 8bit PNGs. One of the CAVE image is saved in 8bpp instead of 16bpp...

// https://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program
bool is_big_endian(void)
{
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}


// Works for 16-bits GRAY PNGs, tested ONLY on a little endian system
// use it at your own risk on big endian system
bool get_png_image_buffer(
    const char* filename,
    std::vector<uint16_t>& framebuffer,
    uint32_t *width, uint32_t *height)
{
    std::FILE * f_png = NULL;
    png_structp png_ptr;
    png_infop info_ptr;
    uint8_t sig[8];
    int bit_depth, color_type;
    std::vector<png_bytep> row_pointers;
#ifndef NDEBUG
    png_uint_32 rowbytes;
#endif // NDEBUG

    f_png = std::fopen(filename, "rb");

    std::fread(sig, 1, 8, f_png);

    if (!png_check_sig(sig, 8)) {
        std::cerr << "Not a PNG!" << std::endl;
        return false;
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr) {
        return false;
    }

    info_ptr = png_create_info_struct(png_ptr);

    if (!info_ptr) {
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        return false;
    }

    png_init_io(png_ptr, f_png);
    png_set_sig_bytes(png_ptr, 8);
    png_read_info(png_ptr, info_ptr);
    png_get_IHDR(png_ptr, info_ptr, width, height, &bit_depth, &color_type, NULL, NULL, NULL);

    // double gamma;
    // if (png_get_gAMA(png_ptr, info_ptr, &gamma)) {
    //     std::cout << "Gamma: " << gamma << std::endl;
    // }

    if (color_type != PNG_COLOR_TYPE_GRAY || bit_depth != 16) {
        std::cerr << "I cannot read this file, must be Gray & 16bit :(" << std::endl;
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        return false;
    }

#ifndef NDEBUG
    // Now we are sure we are in 16bit monochrome
    rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    assert(rowbytes == (*width) * 2);
#endif // NDEBUG

    row_pointers.resize((*height));
    framebuffer.resize((*width) * (*height));

    for (uint32_t y = 0; y < (*height); y++) {
        row_pointers[y] = (png_bytep)&framebuffer[y * (*width)];
    }

    png_read_update_info(png_ptr, info_ptr);
    png_read_image(png_ptr, row_pointers.data());
    png_read_end(png_ptr, NULL);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    // Flip endinaness, remove on big endian processors
    if (!is_big_endian()) {
        for (uint32_t i = 0; i < (*width)*(*height); i++) {
            const int low  = framebuffer[i] & 0xFF;
            const int high = (framebuffer[i] & 0xFF00) >> 8;

            framebuffer[i] = (low << 8) | high;
        }
    } else {
        std::cout << "No byte flipping, executing on a big endian system - NOT TESTED!" << std::endl;
    }

    return true;
}


int main(int argc, char* argv[])
{
    std::string path_first_image, path_output_image;
    Imf::PixelType pixelType = Imf::PixelType::NUM_PIXELTYPES;
    bool write_rgb;

    // Parse arguments
    try {
        TCLAP::CmdLine cmd("Converts a CAVE spectral image to Spectral EXR");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Path to the first PNG image of the CAVE spectral data.", true, "cd_ms_01.png", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Path to the OpenEXR image to write into.", true, "cd_ms.exr", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        std::vector<std::string> allowedPixelTypes;
        allowedPixelTypes.push_back("half");
        allowedPixelTypes.push_back("float");
        allowedPixelTypes.push_back("uint32");

        TCLAP::ValuesConstraint<std::string> allowedPixelTypesVals(allowedPixelTypes);
        TCLAP::ValueArg<std::string> pixelTypeArg("t", "type", "Pixel type to use in the resulting OpenEXR file.", false, "uint32", &allowedPixelTypesVals);

        cmd.add(pixelTypeArg);

        TCLAP::SwitchArg ignoreRgbArg("r", "no_rgb", "Do not write RGB version of the image in the resulting OpenEXR file.", false);

        cmd.add(ignoreRgbArg);

        cmd.parse(argc, argv);

        path_first_image  = inputFileArg.getValue();
        path_output_image = outputFileArg.getValue();

        if (pixelTypeArg.getValue() == "half") {
            pixelType = Imf::PixelType::HALF;
        } else if (pixelTypeArg.getValue() == "float") {
            pixelType = Imf::PixelType::FLOAT;
        } else if (pixelTypeArg.getValue() == "uint32") {
            pixelType = Imf::PixelType::UINT;
        }

        write_rgb = !ignoreRgbArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;

        return 1;
    }

    const std::string path_root = path_first_image.substr(0, path_first_image.size() - 6);

    // Constants for CAVE database
    const int n_bands                 = 31;
    const int start_wavelength_nm     = 400;
    const int wavelength_increment_nm = 10;

    // Read contents of PNGs
    uint32_t width = 0, height = 0;
    bool success = false;

    std::vector<int> wavelengths_nm(n_bands);
    std::vector<std::vector<uint16_t>> original_framebuffers(n_bands);

    for (int band = 0; band < n_bands; band++) {
        uint32_t curr_w, curr_h;

        std::ostringstream ss;
        ss << std::setfill('0') << std::setw(2) << (band + 1);
        std::string current_filename = path_root + ss.str() + ".png";

        success = get_png_image_buffer(current_filename.c_str(), original_framebuffers[band], &curr_w, &curr_h);

        if (!success) {
            std::cerr << "Could not load: " << current_filename << ". Exiting" << std::endl;
            break;
        }

        if (band == 0) {
            width  = curr_w;
            height = curr_h;
        } else {
            if (curr_w != width || curr_h != height) {
                success = false;
                std::cerr << "PNG sizes do not match." << std::endl;
                break;
            }
        }

        wavelengths_nm[band] = start_wavelength_nm + band * wavelength_increment_nm;
    }

    // Convert to OpenEXR
    if (success) {
        // Perform type conversion

        std::vector<float> framebuffers_float;
        std::vector<Imath::half> framebuffers_half;
        std::vector<uint32_t> framebuffers_uint32_t;
        std::vector<float> rgb_images;

        char* data_ptr = NULL;
        size_t type_stride = 0;

        switch (pixelType) {
            case Imf::PixelType::HALF:
                framebuffers_half.resize(width * height * n_bands);

                for (size_t band = 0; band < n_bands; band++) {
                    for (size_t px = 0; px < width * height; px++) {
                        framebuffers_half[px * n_bands + band] = imath_float_to_half(float(original_framebuffers[band][px]) / float(65535));
                    }
                }

                data_ptr = (char*)framebuffers_half.data();
                type_stride = 2;
                break;

            case Imf::PixelType::FLOAT:
                framebuffers_float.resize(width * height * n_bands);
                for (size_t band = 0; band < n_bands; band++) {
                    for (size_t px = 0; px < width * height; px++) {
                        framebuffers_float[px * n_bands + band] = float(original_framebuffers[band][px]) / float(65535);
                    }
                }

                data_ptr = (char*)framebuffers_float.data();
                type_stride = 4;
                break;

            case Imf::PixelType::UINT:
                framebuffers_uint32_t.resize(width * height * n_bands);
                for (size_t band = 0; band < n_bands; band++) {
                    for (size_t px = 0; px < width * height; px++) {
                        framebuffers_uint32_t[px * n_bands + band] = original_framebuffers[band][px];
                    }
                }

                data_ptr = (char*)framebuffers_uint32_t.data();
                type_stride = 4;
                break;

            default:
                std::cerr << "Incorect type selected for writing" << std::endl;
                return 1;
        }

        Imf::Header       exr_header(width, height);
        Imf::ChannelList &exr_channels = exr_header.channels();
        Imf::FrameBuffer  exr_framebuffer;

        SpectrumConverter emissive_converter(true);

        const size_t x_stride = n_bands * type_stride;
        const size_t y_stride = x_stride * width;

        for (size_t band = 0; band < n_bands; band++) {
            const std::string layer_name = "S0." + std::to_string(wavelengths_nm[band]) + "nm";

            exr_channels.insert(layer_name, Imf::Channel(pixelType));

            exr_framebuffer.insert(
                layer_name,
                Imf::Slice(
                    pixelType,
                    &data_ptr[type_stride * band],
                    x_stride, y_stride)
            );
        }

        if (write_rgb) {
            std::vector<float> wavelengths_nm_f(n_bands);

            framebuffers_float.resize(width * height * n_bands);

            for (size_t band = 0; band < n_bands; band++) {
                wavelengths_nm_f[band] = wavelengths_nm[band];

                for (size_t px = 0; px < width * height; px++) {
                    framebuffers_float[px * n_bands + band] = float(original_framebuffers[band][px]) / float(65535);
                }
            }

            emissive_converter.spectralImageToRGB(wavelengths_nm_f, framebuffers_float, width, height, rgb_images);

            const char* RGB[3] = {"R", "G", "B"};

            for (int c = 0; c < 3; c++) {
                exr_channels.insert(RGB[c], Imf::Channel(Imf::FLOAT));

                exr_framebuffer.insert(
                    RGB[c],
                    Imf::Slice(
                        Imf::FLOAT,
                        (char*)&rgb_images[c],
                        3 * sizeof(float),
                        3 * width * sizeof(float))
                );
            }
        }

        Imf::OutputFile exr_out(path_output_image.c_str(), exr_header);
        exr_out.setFrameBuffer(exr_framebuffer);
        exr_out.writePixels(height);
    }

    return 0;
}

/**
 * Copyright 2022 Alban Fichet
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

#include <cstdint>
#include <cmath>
#include <cassert>

#include <png.h>


#include <EXRSpectralImage.h>

// TEMP
// #include <fstream>


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
    png_uint_32 rowbytes;
    std::vector<png_bytep> row_pointers;

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

    // Now we are sure we are in 16bit monochrome
    rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    assert(rowbytes == (*width) * 2);

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
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <path_first_image> <spectral_exr_out>" << std::endl;
    }

    const std::string path_first_image  = argv[1];
    const std::string path_output_image = argv[2];

    const std::string path_root = path_first_image.substr(0, path_first_image.size() - 6);

    const int n_bands = 31;
    const int start_wavelength_nm = 400;
    const int wavelength_increment_nm = 10;

    uint32_t width, height;
    bool success = false;

    std::vector<float> wavelengths_nm(n_bands);
    std::vector<float> spectral_framebuffer;

    for (int i = 0; i < n_bands; i++) {
        uint32_t curr_w, curr_h;
        std::vector<uint16_t> grey_framebuffer;

        std::ostringstream ss;
        ss << std::setfill('0') << std::setw(2) << (i + 1);
        std::string current_filename = path_root + ss.str() + ".png";

        success = get_png_image_buffer(current_filename.c_str(), grey_framebuffer, &curr_w, &curr_h);

        if (!success) {
            std::cerr << "Could not load: " << current_filename << ". Exiting" << std::endl;
            break;
        }

        if (i == 0) {
            width  = curr_w;
            height = curr_h;
            spectral_framebuffer.resize(n_bands * width * height);
        } else {
            if (curr_w != width || curr_h != height) {
                success = false;
                std::cerr << "PNG sizes do not match." << std::endl;
                break;
            }
        }

        wavelengths_nm[i] = start_wavelength_nm + i * wavelength_increment_nm;

        const int stride = n_bands;
        const int start_offset = i;

        for (uint32_t i = 0; i < width * height; i++) {
            assert(start_offset + i * stride < n_bands * width * height);

            spectral_framebuffer[start_offset + i * stride] = 
                (float)grey_framebuffer[i]
                / 65535.f;
        }
    }

    if (success) {
        EXRSpectralImage exr_out(width, height);
        exr_out.appendSpectralFramebuffer(wavelengths_nm, spectral_framebuffer, "S0", PixelType::HALF);

        exr_out.write(path_output_image.c_str());
    }

    return 0;
}
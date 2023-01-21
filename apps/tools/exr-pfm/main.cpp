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
#include <sstream>
#include <cstdint>
#include <cstring>

#include <EXRImage.h>
#include <Util.h>

bool is_big_endian(void)
{
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}


bool write_pfm(
    const char* filename_out,
    const std::vector<float>& framebuffer,
    uint32_t width, uint32_t height)
{
    if (framebuffer.size() != width * height) {
        std::cerr << "Framebuffer dimentions are incorrect" << std::endl;
        return false;
    }

    const char type[2] = {'P', 'f'};
    const char newline = 0x0a;
    const char* little_endian_marker = "-1.0";
    const char* big_endian_marker    = "1.0";

    char width_str[6] = {0};
    char height_str[6] = {0};

    std::sprintf(width_str, "%d ", width);
    std::sprintf(height_str, "%d", height);

    std::FILE* image_out = std::fopen(filename_out, "wb");

    if (!image_out) {
        std::cerr << "Could not open " << filename_out << " for writing." << std::endl;
        return false;
    }

    // Header
    std::fwrite(type, 1, 2, image_out);
    std::fwrite(&newline, 1, 1, image_out);
    std::fwrite(width_str, 1, std::strlen(width_str), image_out);
    std::fwrite(height_str, 1, std::strlen(height_str), image_out);
    std::fwrite(&newline, 1, 1, image_out);

    if (is_big_endian()) {
        std::fwrite(big_endian_marker, 1, std::strlen(big_endian_marker), image_out);
    } else {
        std::fwrite(little_endian_marker, 1, std::strlen(little_endian_marker), image_out);
    }

    std::fwrite(&newline, 1, 1, image_out);

    // Main image
    std::fwrite(framebuffer.data(), sizeof(float), framebuffer.size(), image_out);

    std::fclose(image_out);

    return true;
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "Usage" << std::endl
                  << "-----" << std::endl
                  << argv[0] << " <exr_file_in> <pfm_file_out>" << std::endl;
        return 0;
    }

    const char* filename_in = argv[1];
    const char* filename_out = argv[2];

    const EXRImage image_in(filename_in);

    const uint32_t width  = image_in.width();
    const uint32_t height = image_in.height();

    const std::vector<EXRFramebuffer*>& framebuffers = image_in.getFramebuffersConst();
    const bool has_multiple_framebuffers = framebuffers.size() > 1;

    for (size_t i = 0; i < framebuffers.size(); i++) {
        std::stringstream ss;

        if (has_multiple_framebuffers) {
            std::string base, ext;
            Util::split_extension(filename_out, base, ext);

            ss << base << "_" << framebuffers[i]->getName() << ext;
        } else {
            ss << filename_out;
        }

        if (!write_pfm(ss.str().c_str(), framebuffers[i]->getPixelDataConst(), width, height)) {
            std::cerr << "Failed to write layer " << framebuffers[i]->getName() << std::endl;
        }
    }

    return 0;
}

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

#include <tclap/CmdLine.h>

#include <EXRImage.h>

int main(int argc, char* argv[]) {

    std::string path_image_in, path_image_out;
    Imf::Compression compression;
    bool change_pixel_type = false;
    Imf::PixelType pixel_type;

    try {
        TCLAP::CmdLine cmd("Change the compression used in an OpenEXR image");

        TCLAP::UnlabeledValueArg<std::string> inputFileArg ("Input", "Path to OpenEXR image read from.", true, "in.exr", "path");
        TCLAP::UnlabeledValueArg<std::string> outputFileArg("Output", "Path to the OpenEXR image to write to.", true, "out.exr", "path");

        cmd.add(inputFileArg);
        cmd.add(outputFileArg);

        std::vector<std::string> allowedCompressionConstraint;
        allowedCompressionConstraint.push_back("none");
        allowedCompressionConstraint.push_back("rle");
        allowedCompressionConstraint.push_back("zips");
        allowedCompressionConstraint.push_back("zip");
        allowedCompressionConstraint.push_back("piz");
        allowedCompressionConstraint.push_back("pxr24");
        allowedCompressionConstraint.push_back("b44");
        allowedCompressionConstraint.push_back("b44a");
        allowedCompressionConstraint.push_back("dwaa");
        allowedCompressionConstraint.push_back("dwab");

        TCLAP::ValuesConstraint<std::string> allowedCompressionMethodsVals(allowedCompressionConstraint);
        TCLAP::ValueArg<std::string> compressionMethodArg("m", "method", "Specify the compression method to use for the output OpenEXR image.", true, "zip", &allowedCompressionMethodsVals);

        cmd.add(compressionMethodArg);

        std::vector<std::string> allowedPixelTypeConstraint;
        allowedPixelTypeConstraint.push_back("same");
        allowedPixelTypeConstraint.push_back("float");
        allowedPixelTypeConstraint.push_back("half");
        allowedPixelTypeConstraint.push_back("uint");

        TCLAP::ValuesConstraint<std::string> allowedPixelTypeVals(allowedPixelTypeConstraint);
        TCLAP::ValueArg<std::string> pixelTypeArg("t", "type", "Specify the type to use for framebuffers.", false, "same", &allowedPixelTypeVals);

        cmd.add(pixelTypeArg);

        cmd.parse(argc, argv);

        path_image_in = inputFileArg.getValue();
        path_image_out = outputFileArg.getValue();

        if (compressionMethodArg.getValue() == "none") {
            compression = Imf::Compression::NO_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "rle") {
            compression = Imf::Compression::RLE_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "zips") {
            compression = Imf::Compression::ZIPS_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "zip") {
            compression = Imf::Compression::ZIP_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "pxr24") {
            compression = Imf::Compression::PXR24_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "b44") {
            compression = Imf::Compression::B44_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "b44a") {
            compression = Imf::Compression::B44A_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "dwaa") {
            compression = Imf::Compression::DWAA_COMPRESSION;
        } else if (compressionMethodArg.getValue() == "dwab") {
            compression = Imf::Compression::DWAB_COMPRESSION;
        } else {
            std::cerr << "Error: unknown compression type: " << compressionMethodArg.getValue() << std::endl;
            return 1;
        }

        if (pixelTypeArg.isSet()) {
            if (pixelTypeArg.getValue() == "same") {
                change_pixel_type = false;
            } else if (pixelTypeArg.getValue() == "half") {
                change_pixel_type = true;
                pixel_type = Imf::PixelType::HALF;
            } else if (pixelTypeArg.getValue() == "float") {
                change_pixel_type = true;
                pixel_type = Imf::PixelType::FLOAT;
            } else if (pixelTypeArg.getValue() == "uint") {
                change_pixel_type = true;
                pixel_type = Imf::PixelType::UINT;
            }
        } else {
            change_pixel_type = false;
        }
    }  catch (TCLAP::ArgException  &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        return 1;
    }

    const EXRImage image_in(path_image_in);

    EXRImage image_out(image_in.width(), image_in.height(), compression);

    for (const EXRFramebuffer* fb: image_in.getFramebuffersConst()) {
        Imf::PixelType channel_pixel_type;

        if (change_pixel_type) {
            channel_pixel_type = pixel_type;
        } else {
            channel_pixel_type = fb->getPixelType();
        }

        image_out.appendFramebuffer(fb->getPixelDataConst(), fb->getName(), channel_pixel_type);
    }

    image_out.setAttributesData(image_in.getAttributesData());

    image_out.write(path_image_out);
}

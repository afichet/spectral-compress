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

#include "EXRImage.h"

#include <cstring>
#include <cassert>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include "EXRArrayStream.h"
#include "EXRUtils.h"


EXRFramebuffer::EXRFramebuffer(
    uint32_t width, uint32_t height,
    const char* name,
    Imf::PixelType save_as_type)
    : _pixel_data(width * height)
    , _name(std::strlen(name) + 1)
    , _pixel_type(save_as_type)
{
    if (save_as_type >= Imf::PixelType::NUM_PIXELTYPES) {
        throw std::runtime_error("Invalid type provided");
    }

    memcpy(_name.data(), name, std::strlen(name));
    _name[strlen(name)] = 0;
}


EXRFramebuffer::EXRFramebuffer(
    const std::vector<float>& framebuffer,
    const char* name,
    Imf::PixelType save_as_type)
    : _pixel_data(framebuffer)
    , _name(std::strlen(name) + 1)
    , _pixel_type(save_as_type)
{
    if (save_as_type >= Imf::PixelType::NUM_PIXELTYPES) {
        throw std::runtime_error("Invalid type provided");
    }

    memcpy(_name.data(), name, std::strlen(name));
    _name[strlen(name)] = 0;
}


EXRFramebuffer::~EXRFramebuffer()
{
}

// ----------------------------------------------------------------------------

EXRImage::EXRImage(const char* filename)
{
    load(filename);
}


EXRImage::EXRImage(const std::string& filename)
{
    load(filename.c_str());
}


EXRImage::EXRImage(uint32_t width, uint32_t height, Imf::Compression compression)
    : _width(width)
    , _height(height)
    , _compression(compression)
{
}


EXRImage::~EXRImage()
{
    for (size_t i = 0; i < _framebuffers.size(); i++) {
        delete _framebuffers[i];
    }
}


void EXRImage::appendFramebuffer(
    const std::vector<float>& framebuffer,
    const char* name,
    Imf::PixelType save_as_type)
{
    assert(framebuffer.size() == _width * _height);

    if (save_as_type >= Imf::PixelType::NUM_PIXELTYPES) {
        throw std::runtime_error("Invalid type provided");
    }

    EXRFramebuffer *fb = new EXRFramebuffer(framebuffer, name, save_as_type);

    _framebuffers.push_back(fb);
}


void EXRImage::write(const char* filename) const
{
    Imf::Header       exr_header(_width, _height, 1, Imath::V2f(0, 0), 1, Imf::LineOrder::INCREASING_Y, _compression);
    Imf::ChannelList &exr_channels = exr_header.channels();
    Imf::FrameBuffer  exr_framebuffer;

    // Write attributes if any
    if (_attributes_data.size() > 0) {
        EXRArrayStream attr_stream(_attributes_data);

        do {
            char c;
            std::stringstream attribute_name_stream;
            std::stringstream attribute_type_stream;

            do {
                attr_stream.read(&c, 1);
                attribute_name_stream << c;
            } while(c != 0);

            do {
                attr_stream.read(&c, 1);
                attribute_type_stream << c;
            } while (c != 0);

            const std::string attribute_name = attribute_name_stream.str();
            const std::string attribute_type = attribute_type_stream.str();

            const uint32_t sz = attr_stream.size() - attr_stream.tellg();

            // We saved the string size before the value
            if (std::strcmp(attribute_type.c_str(), "string") == 0) {
                attr_stream.read((char*)&sz, sizeof(uint32_t));
            }

            Imf::Attribute* attr = Imf::Attribute::newAttribute(attribute_type.c_str());
            attr->readValueFrom(attr_stream, sz, 1);

            exr_header.insert(attribute_name, *attr);

            delete attr;
        } while (!(attr_stream.tellg() == attr_stream.size()));
    }

    // Write framebuffers
    std::vector<std::vector<half>> half_fb_data;
    std::vector<std::vector<uint32_t>> int_fb_data;
    std::vector<std::vector<float>> rgb_images;

    for (const EXRFramebuffer* fb: _framebuffers) {
        exr_channels.insert(fb->getName(), Imf::Channel(Imf::FLOAT));

        size_t type_stride;

        char* fb_data = condition_framebuffers(
            fb->getPixelDataConst(),
            fb->getPixelType(),
            type_stride,
            int_fb_data,
            half_fb_data
            );

        const size_t x_stride = type_stride;
        const size_t y_stride = x_stride * _width;

        exr_framebuffer.insert(
            fb->getName(),
            Imf::Slice(
                fb->getPixelType(),
                fb_data,
                x_stride, y_stride)
        );
    }

    Imf::OutputFile exr_out(filename, exr_header);
    exr_out.setFrameBuffer(exr_framebuffer);
    exr_out.writePixels(_height);
}


void EXRImage::write(const std::string& filename) const
{
    write(filename.c_str());
}


void EXRImage::load(const char* filename)
{
    Imf::InputFile exr_in(filename);

    const Imf::Header&  exr_header       = exr_in.header();
    const Imath::Box2i& exr_datawindow   = exr_header.dataWindow();
    const Imf::ChannelList &exr_channels = exr_header.channels();

    _width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    _height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    _compression = exr_header.compression();

    // Read attributes
    EXRArrayStream attr_stream;

    for (Imf::Header::ConstIterator it = exr_header.begin(); it != exr_header.end(); it++) {
        const char* attribute_name = it.name();
        const char* attribute_type = it.attribute().typeName();

        if (std::strcmp(attribute_name, "channels") != 0
         && std::strcmp(attribute_name, "compression") != 0
         && std::strcmp(attribute_name, "lineOrder") != 0) {
            attr_stream.write(attribute_name, std::strlen(attribute_name) + 1);
            attr_stream.write(attribute_type, std::strlen(attribute_type) + 1);

            // For unknown reasons as for now, we need to save the length of the string
            if (std::strcmp(attribute_type, "string") == 0) {
                uint32_t str_sz = ((const Imf::StringAttribute&)it.attribute()).value().size();
                attr_stream.write((char*)&str_sz, sizeof(uint32_t));
            }

            it.attribute().writeValueTo(attr_stream, 1);
        }
    }

    _attributes_data = attr_stream.data();

    // Read framebuffers
    Imf::FrameBuffer exr_framebuffer;

    const size_t x_stride = sizeof(float);
    const size_t y_stride = x_stride * _width;

    for (Imf::ChannelList::ConstIterator channel = exr_channels.begin();
        channel != exr_channels.end();
        channel++) {
        EXRFramebuffer *fb = new EXRFramebuffer(_width, _height, channel.name(), channel.channel().type);

        Imf::Slice slice = Imf::Slice::Make(
            Imf::FLOAT,
            fb->getPixelData().data(),
            exr_header.dataWindow(),
            x_stride, y_stride);

        exr_framebuffer.insert(channel.name(), slice);

        _framebuffers.push_back(fb);
    }

    exr_in.setFrameBuffer(exr_framebuffer);
    exr_in.readPixels(exr_datawindow.min.y, exr_datawindow.max.y);
}

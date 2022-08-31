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

#include "EXRImage.h"

#include <cstring>
#include <cassert>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

// ----------------------------------------------------------------------------

EXRArrayStream::EXRArrayStream()
: Imf::OStream("mem")
, Imf::IStream("mem")
, _pos(0)
{}


EXRArrayStream::EXRArrayStream(const std::vector<uint8_t> data)
: Imf::OStream("mem")
, Imf::IStream("mem")
, _data(data)
, _pos(0)
{}


void EXRArrayStream::write(const char c[/*n*/], int n) {
    const uint64_t remaining_bytes = _data.size() - _pos;

    if (remaining_bytes < n) {
        _data.resize(_data.size() + n - remaining_bytes);
    }

    std::memcpy(&_data[_pos], c, n);

    _pos += n;
}


bool EXRArrayStream::read(char c[/*n*/], int n) {
    const uint64_t remaining_bytes = _data.size() - _pos;

    if (remaining_bytes < n) {
        throw std::exception();
    }

    std::memcpy(c, &_data[_pos], n);

    _pos += n;

    return _pos == _data.size();
}


uint64_t EXRArrayStream::tellp() {
    return _pos;
}


uint64_t EXRArrayStream::tellg() {
    return _pos;
}


void EXRArrayStream::seekp(uint64_t pos) {
    _pos = pos;
}


void EXRArrayStream::seekg(uint64_t pos) {
    _pos = pos;
}


const std::vector<uint8_t>& EXRArrayStream::data() const {
    return _data;
}


size_t EXRArrayStream::size() const {
    return _data.size();
}

// ----------------------------------------------------------------------------

EXRFramebuffer::EXRFramebuffer(
    uint32_t width, uint32_t height,
    const char* name)
    : _pixel_data(width * height)
    , _name(new char[std::strlen(name) + 1])
{
    std::strcpy(_name, name);
    _name[strlen(name)] = 0;
}


EXRFramebuffer::EXRFramebuffer(
    const std::vector<float>& framebuffer,
    const char* name)
    : _pixel_data(framebuffer)
    , _name(new char[std::strlen(name) + 1])
{
    std::strcpy(_name, name);
    _name[strlen(name)] = 0;
}


EXRFramebuffer::~EXRFramebuffer()
{
    delete[] _name;
}

// ----------------------------------------------------------------------------

EXRImage::EXRImage(const char* filename)
{
    Imf::InputFile exr_in(filename);

    const Imf::Header&  exr_header       = exr_in.header();
    const Imath::Box2i& exr_datawindow   = exr_header.dataWindow();
    const Imf::ChannelList &exr_channels = exr_header.channels();

    _width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    _height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    // Read attributes
    EXRArrayStream attr_stream;

    for (Imf::Header::ConstIterator it = exr_header.begin(); it != exr_header.end(); it++) {
        const char* attribute_name = it.name();
        const char* attribute_type = it.attribute().typeName();

        if (std::strcmp(attribute_name, "channels") != 0 
         && std::strcmp(attribute_name, "compression") != 0) {
            attr_stream.write(attribute_name, std::strlen(attribute_name) + 1);
            attr_stream.write(attribute_type, std::strlen(attribute_type) + 1);
            
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
        EXRFramebuffer *fb = new EXRFramebuffer(_width, _height, channel.name());

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


EXRImage::EXRImage(uint32_t width, uint32_t height)
    : _width(width)
    , _height(height)
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
    const char* name)
{
    assert(framebuffer.size() == _width * _height);

    EXRFramebuffer *fb = new EXRFramebuffer(framebuffer, name);

    _framebuffers.push_back(fb);
}


void EXRImage::write(const char* filename) const
{
    Imf::Header       exr_header(_width, _height);
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

            Imf::Attribute* attr = Imf::Attribute::newAttribute(attribute_type.c_str());
            attr->readValueFrom(attr_stream, attr_stream.size() - attr_stream.tellg(), 1);

            if (std::strcmp(attribute_name.c_str(), "channels") != 0) {
                exr_header.insert(attribute_name, *attr);
            }

            delete attr;
        } while (!(attr_stream.tellg() == attr_stream.size()));
    }

    // Write framebuffers
    const size_t x_stride = sizeof(float);
    const size_t y_stride = x_stride * _width;

    for (EXRFramebuffer* fb: _framebuffers) {
        exr_channels.insert(fb->getName(), Imf::Channel(Imf::FLOAT));
        
        exr_framebuffer.insert(
            fb->getName(), 
            Imf::Slice(
                Imf::FLOAT, 
                (char*)fb->getPixelData().data(), 
                x_stride, y_stride)
        );
    }

    Imf::OutputFile exr_out(filename, exr_header);
    exr_out.setFrameBuffer(exr_framebuffer);
    exr_out.writePixels(_height);
}
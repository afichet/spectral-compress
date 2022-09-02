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

#include "EXRSpectralImage.h"

#include <cstring>
#include <cassert>
#include <map>
#include <set>
#include <regex>
#include <algorithm>
#include <exception>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include "EXRArrayStream.h"

#include "SpectrumConverter.h"

void vec_float_to_half(
    const std::vector<float>& framebuffer_in,
    std::vector<half> & framebuffer_out)
{
    framebuffer_out.resize(framebuffer_in.size());

    for (size_t i = 0; i < framebuffer_in.size(); i++) {
        float v = framebuffer_in[i];
        if (v > 65535.f) v = 65535.f;

        framebuffer_out[i] = (half)v;
    }
}


void vec_float_to_int32(
    const std::vector<float>& framebuffer_in,
    std::vector<uint32_t>& framebuffer_out)
{
    framebuffer_out.resize(framebuffer_in.size());

    for (size_t i = 0; i < framebuffer_in.size(); i++) {
        float v = framebuffer_in[i];
        if (v > 1.f) v = 1.f;

        framebuffer_out[i] = (uint32_t)(v * 65535.f);
    }
}


PixelType from_Imf(Imf::PixelType pixel_type)
{
    switch (pixel_type) {
        case Imf::PixelType::UINT:  return PixelType::UINT;
        case Imf::PixelType::HALF:  return PixelType::HALF;
        case Imf::PixelType::FLOAT: return PixelType::FLOAT;
        default: 
            throw std::runtime_error("Pixel type is not handled!");
    }
}


Imf::PixelType to_Imf(PixelType pixel_type)
{
    switch (pixel_type) {
        case PixelType::UINT:  return Imf::PixelType::UINT;
        case PixelType::HALF:  return Imf::PixelType::HALF;
        case PixelType::FLOAT: return Imf::PixelType::FLOAT;
        default: 
            throw std::runtime_error("Pixel type is not handled!");
    }
}

char* condition_framebuffers(
    const std::vector<float>& framebuffer_in,
    Imf::PixelType write_pixel_type,
    size_t& type_stride,
    std::vector<std::vector<uint32_t>>& framebuffers_int,
    std::vector<std::vector<half>>& framebuffers_half) 
{
    char* fb_data;

    switch (write_pixel_type)
    {
        case Imf::PixelType::UINT:
            type_stride = sizeof(uint32_t);
            framebuffers_int.resize(framebuffers_int.size() + 1);
            vec_float_to_int32(framebuffer_in, framebuffers_int.back());
            fb_data = (char*)framebuffers_int.back().data();
            break;
        
        case Imf::PixelType::HALF:
            type_stride = sizeof(half);
            framebuffers_half.resize(framebuffers_half.size() + 1);
            vec_float_to_half(framebuffer_in, framebuffers_half.back());
            fb_data = (char*)framebuffers_half.back().data();
            break;

        case Imf::PixelType::FLOAT:
            type_stride = sizeof(float);
            fb_data = (char*)framebuffer_in.data();
            break;

        default:
            throw std::runtime_error("Unknonw pixel_type");
    }

    return fb_data;
}



EXRSpectralImage::EXRSpectralImage(const char* filename)
{
    Imf::InputFile exr_in(filename);

    const Imf::Header&  exr_header       = exr_in.header();
    const Imath::Box2i& exr_datawindow   = exr_header.dataWindow();
    const Imf::ChannelList &exr_channels = exr_header.channels();

    _width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    _height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    // =======================================================================
    // Read attributes
    // =======================================================================

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

    // ========================================================================
    // Read framebuffers
    // =======================================================================
    
    std::map<std::string, std::vector<std::pair<std::string, float>>> spectral_channels;
    
    std::set<std::string> extra_channels;
    std::set<std::string> ignored_channels;

    std::map<std::string, PixelType> channel_type;

    Imf::FrameBuffer exr_framebuffer;

    // ------------------------------------------------------------------------
    // Determine channels' position
    // ------------------------------------------------------------------------

    const std::regex expr(
        "^(.*)((S([0-3]))|T)\\.(\\d*,?\\d*([Ee][+-]?\\d+)?)(Y|Z|E|P|T|G|M|k|h|"
        "da|d|c|m|u|n|p)?(m|Hz)$");

    for (Imf::ChannelList::ConstIterator channel = exr_channels.begin();
            channel != exr_channels.end();
            channel++) {
        const std::string name          = channel.name();

        // seems useless but allows to have a clean .h without any OpenEXR imclude
        const PixelType p_type = from_Imf(channel.channel().type);

        std::smatch matches;
        const bool matched = std::regex_search(name, matches, expr);

        if (matched) {
            const std::string root   = matches[1].str();
            const std::string prefix = root + matches[2].str();
            
            if (spectral_channels[prefix].size() == 0) {
                ignored_channels.insert(root + "R");
                ignored_channels.insert(root + "G");
                ignored_channels.insert(root + "B");                
            }

            const double value_nm = to_nm(
                matches[5].str(),
                matches[7].str(),
                matches[8].str()
            );

            spectral_channels[prefix].push_back(std::make_pair(channel.name(), value_nm));

            channel_type[prefix] = p_type;
        } else {
            extra_channels.insert(name);
            channel_type[name] = p_type;
        }
    }

    // Filter out ignored channels
    for (const std::string& s: ignored_channels) {
        extra_channels.erase(s);
    }

    // ------------------------------------------------------------------------
    // Read framebuffers
    // ------------------------------------------------------------------------

    for (const auto& n: spectral_channels) {
        const std::string& root_name = n.first;
        const std::vector<std::pair<std::string, float>>& wavelengths = n.second;

        SpectralFramebuffer* fb = new SpectralFramebuffer;
        fb->root_name = root_name;
        fb->image_data.resize(_width * _height * wavelengths.size());
        fb->wavelengths_nm.reserve(wavelengths.size());
        fb->pixel_type = channel_type[root_name];

        const size_t x_stride = sizeof(float) * wavelengths.size();
        const size_t y_stride = x_stride * _width;

        for (size_t wl_idx = 0; wl_idx < wavelengths.size(); wl_idx++) {
            const std::string& layer_name    = n.second[wl_idx].first;
            const float        wavelength_nm = n.second[wl_idx].second;

            Imf::Slice slice = Imf::Slice::Make(
                Imf::FLOAT,
                &(fb->image_data[wl_idx]),
                exr_header.dataWindow(),
                x_stride, y_stride);
            
            exr_framebuffer.insert(layer_name, slice);

            fb->wavelengths_nm.push_back(wavelength_nm);
        }

        _spectral_framebuffers.push_back(fb);
    }

    for (const auto& name: extra_channels) {
        GreyFramebuffer* fb = new GreyFramebuffer;
        
        fb->layer_name = name;
        fb->image_data.resize(_width * _height);
        fb->pixel_type = channel_type[name];

        Imf::Slice slice = Imf::Slice::Make(
            Imf::FLOAT,
            fb->image_data.data(),
            exr_header.dataWindow());
            
        exr_framebuffer.insert(name, slice);

        _extra_framebuffers.push_back(fb);
    }

    exr_in.setFrameBuffer(exr_framebuffer);
    exr_in.readPixels(exr_datawindow.min.y, exr_datawindow.max.y);
}


EXRSpectralImage::EXRSpectralImage(uint32_t width, uint32_t height)
    : _width(width)
    , _height(height)
{   
}


EXRSpectralImage::~EXRSpectralImage()
{
    for (SpectralFramebuffer* framebuffer: _spectral_framebuffers) {
        delete framebuffer;
    }

    for (GreyFramebuffer* framebuffer: _extra_framebuffers) {
        delete framebuffer;
    }
}


void EXRSpectralImage::appendSpectralFramebuffer(
    const std::vector<float> &wavelengths_nm,
    const std::vector<float> &framebuffer,
    const std::string& prefix,
    const PixelType save_as_type)
{
    appendSpectralFramebuffer(
        wavelengths_nm,
        framebuffer,
        prefix.c_str(),
        save_as_type
    );
}


void EXRSpectralImage::appendSpectralFramebuffer(
    const std::vector<float> &wavelengths_nm,
    const std::vector<float> &framebuffer,
    const char* prefix,
    const PixelType save_as_type)
{
    assert(framebuffer.size() == _width * _height * wavelengths_nm.size());

    if (save_as_type >= N_PIXEL_TYPES) {
        throw std::runtime_error("Invalid type provided");
    }

    SpectralFramebuffer* fb = new SpectralFramebuffer;
    
    fb->root_name      = prefix;
    fb->wavelengths_nm = wavelengths_nm;
    fb->image_data     = framebuffer;
    fb->pixel_type     = save_as_type;

    _spectral_framebuffers.push_back(fb);
}


void EXRSpectralImage::appendExtraFramebuffer(
    const std::vector<float>& framebuffer,
    const std::string& name,
    const PixelType save_as_type)
{
    appendExtraFramebuffer(
        framebuffer,
        name.c_str(),
        save_as_type
    );
}


void EXRSpectralImage::appendExtraFramebuffer(
    const std::vector<float>& framebuffer,
    const char* name,
    const PixelType save_as_type)
{
    assert(framebuffer.size() == _width * _height);

    if (save_as_type >= N_PIXEL_TYPES) {
        throw std::runtime_error("Invalid type provided");
    }

    GreyFramebuffer* fb = new GreyFramebuffer;

    fb->layer_name = name;
    fb->image_data = framebuffer;
    fb->pixel_type = save_as_type;

    _extra_framebuffers.push_back(fb);
}


inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}


void EXRSpectralImage::write(const char* filename) const
{
    Imf::Header       exr_header(_width, _height);
    Imf::ChannelList &exr_channels = exr_header.channels();
    Imf::FrameBuffer  exr_framebuffer;

    SpectrumConverter reflective_converter(false);
    SpectrumConverter emissive_converter(true);

    // ------------------------------------------------------------------------
    // Write attributes if any
    // ------------------------------------------------------------------------

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

    // ------------------------------------------------------------------------
    // Write framebuffers
    // ------------------------------------------------------------------------

    // Need to keep buffers alive when converting from one type to another
    std::vector<std::vector<half>> half_fb_data;
    std::vector<std::vector<uint32_t>> int_fb_data;
    std::vector<std::vector<float>> rgb_images;

    // Spectral groups
    for (const SpectralFramebuffer* framebuffer: _spectral_framebuffers) {
        const size_t n_bands            = framebuffer->wavelengths_nm.size();
        const Imf::PixelType pixel_type = to_Imf(framebuffer->pixel_type);

        std::string root_name = framebuffer->root_name;

        size_t type_stride;

        char* fb_data = condition_framebuffers(
            framebuffer->image_data,
            pixel_type,
            type_stride,
            int_fb_data,
            half_fb_data
            );

        const size_t x_stride = n_bands * type_stride;
        const size_t y_stride = x_stride * _width;

        for (size_t i = 0; i < n_bands; i++) {
            std::string wl_nm = std::to_string(framebuffer->wavelengths_nm[i]) + "nm";
            std::replace(wl_nm.begin(), wl_nm.end(), '.', ',');
            const std::string layer_name = root_name + "." + wl_nm;

            exr_channels.insert(layer_name, Imf::Channel(pixel_type));
            
            exr_framebuffer.insert(
                layer_name, 
                Imf::Slice(
                    pixel_type,
                    &fb_data[type_stride * i],
                    x_stride, y_stride)
            );
        }

        // Rebuild RGB layers from spectral hierachy
        // TODO: consider the case where the hierachy both has T and S0 children
        bool hasRGB = false;

        if (ends_with(root_name, "S0")) {
            rgb_images.resize(rgb_images.size() + 1);

            emissive_converter.spectralImageToRGB(
                framebuffer->wavelengths_nm,
                framebuffer->image_data,
                _width, _height,
                rgb_images.back()
            );

            hasRGB = true;
        } else if (ends_with(root_name, "T")) {
            rgb_images.resize(rgb_images.size() + 1);

            reflective_converter.spectralImageToRGB(
                framebuffer->wavelengths_nm,
                framebuffer->image_data,
                _width, _height,
                rgb_images.back()
            );

            hasRGB = true;
        }

        if (hasRGB) {
            char* fb_data = condition_framebuffers(
                rgb_images.back(),
                pixel_type,
                type_stride,
                int_fb_data,
                half_fb_data
                );
                
            const char* RGB[3] = {"R", "G", "B"};

            for (int c = 0; c < 3; c++) {
                exr_channels.insert(RGB[c], Imf::Channel(pixel_type));

                exr_framebuffer.insert(
                    RGB[c], 
                    Imf::Slice(
                        pixel_type,
                        &fb_data[c * type_stride],
                        3 * type_stride,
                        3 * _width * type_stride)
                );
            }
        }
    }

    // Gray framebuffers
    for (const GreyFramebuffer* framebuffer: _extra_framebuffers) {
        const Imf::PixelType pixel_type = to_Imf(framebuffer->pixel_type);
        
        size_t type_stride;

        char* fb_data = condition_framebuffers(
            framebuffer->image_data,
            pixel_type,
            type_stride,
            int_fb_data,
            half_fb_data
        );
        
        const size_t x_stride = type_stride;
        const size_t y_stride = x_stride * _width;

        exr_channels.insert(framebuffer->layer_name, Imf::Channel(pixel_type));
        
        exr_framebuffer.insert(
            framebuffer->layer_name,
            Imf::Slice(
                pixel_type,
                fb_data,
                x_stride, y_stride)
        );
    }

    Imf::OutputFile exr_out(filename, exr_header);
    exr_out.setFrameBuffer(exr_framebuffer);
    exr_out.writePixels(_height);
}


double EXRSpectralImage::to_nm(
    const std::string& value,
    const std::string& prefix,
    const std::string& units) 
{
    std::string centralValueStr(value);
    std::replace(
        centralValueStr.begin(),
        centralValueStr.end(),
        ',',
        '.');

    const double v = std::stod(centralValueStr);

    if (prefix == "n" && units == "m") return v;

    double wavelength_nm = v;

    const std::map<std::string, double> unit_prefix = {
        {"Y", 1e24},
        {"Z", 1e21},
        {"E", 1e18},
        {"P", 1e15},
        {"T", 1e12},
        {"G", 1e9},
        {"M", 1e6},
        {"k", 1e3},
        {"h", 1e2},
        {"da", 1e1},
        {"d", 1e-1},
        {"c", 1e-2},
        {"m", 1e-3},
        {"u", 1e-6},
        {"n", 1e-9},
        {"p", 1e-12}};

    // Apply multiplier
    if (prefix.size() > 0) {
        wavelength_nm *= unit_prefix.at(prefix);
    }

    // Apply units
    if (units == "Hz") {
        wavelength_nm = 299792458. / wavelength_nm * 1e9;
    } else if (units == "m") {
        wavelength_nm = wavelength_nm * 1e9;
    } else {
        // Unknown unit
        // Something went wrong with the parsing. This shall not occur.
        throw std::out_of_range("Unknown unit");
    }

    return wavelength_nm;
}

// ----------------------------------------------------------------------------
// Getter setters
// ----------------------------------------------------------------------------

uint32_t EXRSpectralImage::width()  const
{
    return _width;
}


uint32_t EXRSpectralImage::height() const 
{ 
    return _height;
}


std::vector<SpectralFramebuffer*>& EXRSpectralImage::getSpectralFramebuffers() 
{
    return _spectral_framebuffers;
}


std::vector<GreyFramebuffer*>& EXRSpectralImage::getExtraFramebuffers() 
{
    return _extra_framebuffers;
}


void EXRSpectralImage::setAttributesData(const std::vector<uint8_t>& data) 
{ 
    _attributes_data = data; 
}


const std::vector<uint8_t>& EXRSpectralImage::getAttributesData() const 
{ 
    return _attributes_data; 
}


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
#include <regex>
#include <string>
#include <map>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include <OpenEXR/ImfAttribute.h>
#include <OpenEXR/ImfIO.h>

double to_nm(
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


struct SpectralFramebuffer 
{
    std::string root_name;
    std::vector<float> wavelengths_nm;
    std::vector<float> image_data;
};


struct GreyFramebuffer
{
    std::string layer_name;
    std::vector<float> image_data;
};


class ArrayStream : public Imf::OStream, public Imf::IStream {
public:
    ArrayStream()
    : OStream("mem")
    , IStream("mem")
    , _pos(0)
    {}

    
    ArrayStream(const std::vector<uint8_t> data)
    : OStream("mem")
    , IStream("mem")
    , _content(data)
    , _pos(0)
    {}


    virtual void write(const char c[/*n*/], int n) {
        const uint64_t remaining_bytes = _content.size() - _pos;

        if (remaining_bytes < n) {
            _content.resize(_content.size() + n - remaining_bytes);
        }

        std::memcpy(&_content[_pos], c, n);

        _pos += n;
    }

    
    virtual bool read(char c[/*n*/], int n) {
        const uint64_t remaining_bytes = _content.size() - _pos;

        if (remaining_bytes < n) {
            throw std::exception();
        }

        std::memcpy(c, &_content[_pos], n);

        _pos += n;

        return _pos == _content.size();
    }

    
    virtual uint64_t tellp () {
        return _pos;
    }

    
    virtual uint64_t tellg() {
        return _pos;
    }

    
    virtual void seekp(uint64_t pos) {
        _pos = pos;
    }

    
    virtual void seekg(uint64_t pos) {
        _pos = pos;
    }

private:
    std::vector<uint8_t> _content;
    uint64_t _pos;
};


int main(int argc, char *argv[])
{
    Imf::InputFile      exr_in(argv[1]);
    const Imf::Header&  exr_header     = exr_in.header();
    const Imath::Box2i& exr_datawindow = exr_header.dataWindow();

    uint32_t width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    uint32_t height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    ArrayStream stream;

    for (Imf::Header::ConstIterator it = exr_header.begin(); it != exr_header.end(); it++) {
        std::cout << it.name() << " " << std::endl;//<< it.attribute() << std::endl;

        it.attribute().writeValueTo(stream, 1);
    }

    // ---------------------------------------------------------------------
    // Determine channels' position
    // ---------------------------------------------------------------------

    // std::array<std::vector<std::pair<float, std::string>>, 4> wavelengths_nm_S;
    // std::vector<std::pair<float, std::string>> wavelengths_nm_reflective;

    const std::regex expr(
        "^(.*)((S([0-3]))|T)\\.(\\d*,?\\d*([Ee][+-]?\\d+)?)(Y|Z|E|P|T|G|M|k|h|"
        "da|d|c|m|u|n|p)?(m|Hz)$");
    
    // map for keys wl spectral channels
    // map for other framebuffers

    std::map<std::string, std::vector<std::pair<std::string, float>>> spectral_channels;
    std::set<std::string> extra_channels;

    std::set<std::string> ignored_channels;
        
    const Imf::ChannelList &exr_channels = exr_header.channels();

    for (Imf::ChannelList::ConstIterator channel = exr_channels.begin();
            channel != exr_channels.end();
            channel++) {
            // Check if the channel is a spectral one

        std::smatch matches;
        const std::string name = channel.name();
        const bool matched = std::regex_search(name, matches, expr);

        if (matched) {
            std::string root   = matches[1].str();
            std::string prefix = root + matches[2].str();
            
            std::cout << "matched: " << prefix << std::endl;

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
        } else {
            extra_channels.insert(channel.name());
        }
    }

    // filter out ignored channels
    for (const std::string& s: ignored_channels) {
        extra_channels.erase(s);
    }
    
    std::cout << std::endl;

    // 
    std::cout << "Spectral channels:" << std::endl;
    std::cout << "------------------" << std::endl;
    
    for (const auto& n : spectral_channels) {
        std::cout << n.first << "\t[";
        for (const auto& wl: n.second) {
            std::cout << wl.first << ",";
        }
        std::cout << "]" << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Extra channels:" << std::endl;
    std::cout << "---------------" << std::endl;

    for (const auto& n: extra_channels) {
        std::cout << n << std::endl;
    }

    std::cout << spectral_channels.size() << " " << extra_channels.size() << std::endl;

    // Now we can read framebuffers
    std::vector<SpectralFramebuffer*> spectral_framebuffers;
    std::vector<GreyFramebuffer*> extra_framebuffers;

    Imf::FrameBuffer exrFrameBuffer;

    for (const auto& n: spectral_channels) {
        const std::string& root_name = n.first;
        const std::vector<std::pair<std::string, float>>& wavelengths = n.second;

        SpectralFramebuffer* fb = new SpectralFramebuffer;
        fb->root_name = root_name;
        fb->image_data.resize(width * height * wavelengths.size());
        fb->wavelengths_nm.reserve(wavelengths.size());

        const size_t x_stride = sizeof(float) * wavelengths.size();
        const size_t y_stride = x_stride * width;

        for (size_t wl_idx = 0; wl_idx < wavelengths.size(); wl_idx++) {
            const std::string& layer_name    = n.second[wl_idx].first;
            const float        wavelength_nm = n.second[wl_idx].second;

            Imf::Slice slice = Imf::Slice::Make(
                Imf::FLOAT,
                &(fb->image_data[wl_idx]),
                exr_header.dataWindow(),
                x_stride, y_stride);
            
            exrFrameBuffer.insert(layer_name, slice);
        }

        spectral_framebuffers.push_back(fb);
    }

    for (const auto& name: extra_channels) {
        GreyFramebuffer* fb = new GreyFramebuffer;
        
        fb->layer_name = name;

        Imf::Slice slice = Imf::Slice::Make(
            Imf::FLOAT,
            fb->image_data.data(),
            exr_header.dataWindow());
            
        exrFrameBuffer.insert(name, slice);        
    }

    exr_in.setFrameBuffer(exrFrameBuffer);
    exr_in.readPixels(exr_datawindow.min.y, exr_datawindow.max.y);

    return 0;
}
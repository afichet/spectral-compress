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
#include <sstream>
#include <limits>
#include <regex>
#include <string>
#include <map>
#include <set>
#include <algorithm>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include <JXLImage.h>
#include <EXRImage.h>

#include <moments.h>
#include <moments_image.h>


void check_for_invalid_vals(const std::vector<float>& buffer) {
    for (size_t i = 0; i < buffer.size(); i++) {
        const float v = buffer[i];
        if (std::isinf(v) || std::isnan(v)) {
            std::cout << "Invalid moment value!" << std::endl;
        }
    }
}


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


void get_buffers_from_exr(
    Imf::InputFile& exr_in,
    std::vector<SpectralFramebuffer*>& spectral_framebuffers,
    std::vector<GreyFramebuffer*>& extra_framebuffers)
{
    const Imf::Header&  exr_header       = exr_in.header();
    const Imath::Box2i& exr_datawindow   = exr_header.dataWindow();
    const Imf::ChannelList &exr_channels = exr_header.channels();

    const uint32_t width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    const uint32_t height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    std::map<std::string, std::vector<std::pair<std::string, float>>> spectral_channels;
    std::set<std::string> extra_channels;
    std::set<std::string> ignored_channels;

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
        std::smatch matches;
        const std::string name = channel.name();
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
        } else {
            extra_channels.insert(channel.name());
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
            
            exr_framebuffer.insert(layer_name, slice);

            fb->wavelengths_nm.push_back(wavelength_nm);
        }

        spectral_framebuffers.push_back(fb);
    }

    for (const auto& name: extra_channels) {
        GreyFramebuffer* fb = new GreyFramebuffer;
        
        fb->layer_name = name;
        fb->image_data.resize(width * height);

        Imf::Slice slice = Imf::Slice::Make(
            Imf::FLOAT,
            fb->image_data.data(),
            exr_header.dataWindow());
            
        exr_framebuffer.insert(name, slice);

        extra_framebuffers.push_back(fb);
    }

    exr_in.setFrameBuffer(exr_framebuffer);
    exr_in.readPixels(exr_datawindow.min.y, exr_datawindow.max.y);
}


void compress_spectral_framebuffer(
    const SpectralFramebuffer* framebuffer,
    std::vector<std::vector<float>>& compressed_moments,
    std::vector<float>& mins,
    std::vector<float>& maxs)
{
    std::vector<float> phases;
    std::vector<float> moments_image;
    std::vector<float> compressed_moments_image;
    
    const uint32_t n_moments = framebuffer->wavelengths_nm.size();
    const uint32_t n_pixels = framebuffer->image_data.size() / framebuffer->wavelengths_nm.size();

    // TODO: this shall be revisited
    std::vector<float> spectral_framebuffer(framebuffer->image_data.size());

    #pragma omp parallel for
    for (size_t i = 0; i < spectral_framebuffer.size(); i++) {
        float v = framebuffer->image_data[i];

        if (std::isinf(v) || std::isnan(v) || v < 1e-8) {
            v = 1e-8;
        }

        spectral_framebuffer[i] = v;
    }

    wavelengths_to_phases(framebuffer->wavelengths_nm, phases);
    
    compute_moments_image(
        phases, 
        spectral_framebuffer,
        n_pixels, 1, 
        n_moments, 
        moments_image
    );

    compress_moments_image(
        moments_image,
        n_pixels, 1,
        n_moments,
        compressed_moments_image
    );

    compressed_moments.resize(n_moments);

    for (size_t m = 0; m < compressed_moments.size(); m++) {
        compressed_moments[m].resize(n_pixels);
    }

    // DC component does not need further modification
    for (size_t i = 0; i < n_pixels; i++) {
        compressed_moments[0][i] = compressed_moments_image[n_moments * i];
    }

    // Rescale AC components in [0..1]
    for (size_t m = 1; m < n_moments; m++) {
        // Get min / max
        float v_min = compressed_moments_image[m];
        float v_max = compressed_moments_image[m];

        for (size_t i = 0; i < n_pixels; i++) {
            v_min = std::min(v_min, compressed_moments_image[n_moments * i + m]);
            v_max = std::max(v_max, compressed_moments_image[n_moments * i + m]);
        }

        mins.push_back(v_min);
        maxs.push_back(v_max);

        // Now rescale moments
        for (size_t i = 0; i < n_pixels; i++) {
            const float v = compressed_moments_image[n_moments * i + m];

            compressed_moments[m][i] = (v - v_min) / (v_max - v_min);
        }
    }
}


int main(int argc, char *argv[]) 
{
    // TODO: generate a default output name and check if exists,
    // this will allow drag and drop of an EXR over the executable
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <exr_in> <jxl_out>" << std::endl;

        exit(0);
    }

    const bool saveRGB = false;

    const char* filename_in  = argv[1];
    const char* filename_out = argv[2];

    // TODO: move handling of layer types & attributes in EXRImage
    Imf::InputFile exr_in(filename_in);
    Imf::Header exr_header = exr_in.header();

    const Imath::Box2i& dw = exr_in.header().dataWindow();
    const uint32_t width   = dw.max.x - dw.min.x + 1;
    const uint32_t height  = dw.max.y - dw.min.y + 1;

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

    std::vector<SpectralFramebuffer*> spectral_framebuffers;
    std::vector<GreyFramebuffer*> extra_framebuffers;

    get_buffers_from_exr(exr_in, spectral_framebuffers, extra_framebuffers);

    JXLImage jxl_out(width, height);
    SGEGBox box;

    box.exr_attributes = attr_stream.data();
    
    for (const SpectralFramebuffer* fb: spectral_framebuffers) {
        SGEGSpectralGroup sg;

        sg.root_name.resize(fb->root_name.size() + 1);
        std::memcpy(sg.root_name.data(), fb->root_name.c_str(), sg.root_name.size() * sizeof(char));
        sg.wavelengths = fb->wavelengths_nm;

        std::vector<std::vector<float>> compressed_moments;

        compress_spectral_framebuffer(fb, compressed_moments, sg.mins, sg.maxs);

        // Now we can save to JPEG XL
        for (size_t m = 0; m < compressed_moments.size(); m++) {
            float n_bits;
            float n_exponent_bits;

            if (m == 0) {
                n_bits = 32;
                n_exponent_bits = 8;
            } else {
                n_bits = 8;
                n_exponent_bits = 0;
            }

            const size_t idx = jxl_out.appendFramebuffer(
                compressed_moments[m], 
                1, 
                n_bits, 
                n_exponent_bits, 
                1, 
                fb->root_name.c_str());            
            
            sg.layer_indices.push_back(idx);
        }

        box.spectral_groups.push_back(sg);
    }

    for (const GreyFramebuffer* fb: extra_framebuffers) {
        SGEGGrayGroup gg;
        
        gg.layer_name.resize(fb->layer_name.size() + 1);
        std::memcpy(gg.layer_name.data(), fb->layer_name.c_str(), gg.layer_name.size() * sizeof(char));

        // TODO get the original data type from the EXR file
        gg.layer_index = jxl_out.appendFramebuffer(fb->image_data, 1, 32, 8, 1, fb->layer_name.c_str());

        box.gray_groups.push_back(gg);
    }

    // box.print();

    jxl_out.setBox(box);
    jxl_out.write(filename_out);

    for (SpectralFramebuffer* fb: spectral_framebuffers) {
        delete fb;
    }

    for (GreyFramebuffer* fb: extra_framebuffers) {
        delete fb;
    }

    return 0;
}
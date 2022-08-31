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

#include "SGEGBox.h"

#include <cstring>
#include <cstdint>

// ----------------------------------------------------------------------------

size_t SGEGSpectralGroup::getRaw(uint8_t data[]) const
{
    size_t offset = 0, write_size = 0;
    uint32_t data_len = 0;

    data_len = root_name.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(char);
    std::memcpy(data + offset, root_name.data(), write_size);
    offset += write_size;

    // layer_indices
    data_len = layer_indices.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(uint32_t);
    std::memcpy(data + offset, layer_indices.data(), write_size);
    offset += write_size;

    // wavelengths
    data_len = wavelengths.size();
    
    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(float);
    std::memcpy(data + offset, wavelengths.data(), write_size);
    offset += write_size;

    // mins / maxs
    data_len = mins.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(float);
    std::memcpy(data + offset, mins.data(), write_size);
    offset += write_size;
    std::memcpy(data + offset, maxs.data(), write_size);
    offset += write_size;

    return offset;
}


size_t SGEGSpectralGroup::fromRaw(const uint8_t data[])
{
    size_t offset = 0, read_size = 0;
    uint32_t data_len = 0;

    // root_name
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    root_name.resize(data_len);

    read_size = data_len * sizeof(char);
    std::memcpy(root_name.data(), data + offset, read_size);
    offset += read_size;

    // layer_indices
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    layer_indices.resize(data_len);

    read_size = data_len * sizeof(uint32_t);
    std::memcpy(layer_indices.data(), data + offset, read_size);
    offset += read_size;

    // wavelengths
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    wavelengths.resize(data_len);

    read_size = data_len * sizeof(float);
    std::memcpy(wavelengths.data(), data + offset, read_size);
    offset += read_size;

    // mins / maxs
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    mins.resize(data_len);
    maxs.resize(data_len);

    read_size = data_len * sizeof(float);
    std::memcpy(mins.data(), data + offset, read_size);
    offset += read_size;
    std::memcpy(maxs.data(), data + offset, read_size);
    offset += read_size;

    return offset;
}


size_t SGEGSpectralGroup::size() const
{
    return
        sizeof(uint32_t) // Array size
        + root_name.size() * sizeof(char)
        + sizeof(uint32_t) // Array size
        + layer_indices.size() * sizeof(uint32_t)
        + sizeof(uint32_t) // Array size
        + wavelengths.size() * sizeof(float)
        + sizeof(uint32_t) // Arrays size
        + mins.size() * sizeof(float)
        + maxs.size() * sizeof(float);
}


SGEGSpectralGroup& SGEGSpectralGroup::operator=(const SGEGSpectralGroup& other) {
    if (this != &other) {
        root_name     = other.root_name;
        layer_indices = other.layer_indices;
        wavelengths   = other.wavelengths;
        mins          = other.mins;
        maxs          = other.maxs;
    }

    return *this;
}

// ----------------------------------------------------------------------------

size_t SGEGGrayGroup::getRaw(uint8_t data[]) const
{
    size_t offset = 0, write_size = 0;
    uint32_t data_len = 0;

    // layer_name
    data_len = layer_name.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(char);
    std::memcpy(data + offset, layer_name.data(), write_size);
    offset += write_size;

    // layer_index
    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &layer_index, write_size);
    offset += write_size;

    return offset;
}


size_t SGEGGrayGroup::fromRaw(const uint8_t data[])
{
    size_t offset = 0, read_size = 0;
    uint32_t data_len = 0;

    // layer_name
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    layer_name.resize(data_len);

    read_size = data_len * sizeof(char);
    std::memcpy(layer_name.data(), data + offset, read_size);
    offset += read_size;

    // layer_index
    read_size = sizeof(uint32_t);
    std::memcpy(&layer_index, data + offset, read_size);
    offset += read_size;

    return offset;
}


size_t SGEGGrayGroup::size() const
{
    return
        sizeof(uint32_t) // Array size
        + layer_name.size() * sizeof(char)
        + sizeof(uint32_t);
}


SGEGGrayGroup& SGEGGrayGroup::operator=(const SGEGGrayGroup& other) {
    if (this != &other) {
        layer_name  = other.layer_name;
        layer_index = other.layer_index;
    }

    return *this;
}

// ----------------------------------------------------------------------------

SGEGBox::SGEGBox() 
    : revision(1)
{}


SGEGBox::SGEGBox(const std::vector<uint8_t> &data)
{
    uint32_t data_len = 0;
    size_t offset = 0, read_size = 0;

    // revision
    read_size = sizeof(uint32_t);
    std::memcpy(&revision, data.data() + offset, read_size);
    offset += read_size;

    // spectral_groups
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data.data() + offset, read_size);
    offset += read_size;

    spectral_groups.resize(data_len);

    for (size_t i = 0; i < spectral_groups.size(); i++) {
        offset += spectral_groups[i].fromRaw(data.data() + offset);
    } 

    // gray_groups
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data.data() + offset, read_size);
    offset += read_size;

    gray_groups.resize(data_len);

    for (size_t i = 0; i < gray_groups.size(); i++) {
        offset += gray_groups[i].fromRaw(data.data() + offset);
    }

    // exr attributes
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data.data() + offset, read_size);
    offset += read_size;

    exr_attributes.resize(data_len);

    read_size = data_len * sizeof(uint8_t);
    std::memcpy(exr_attributes.data(), data.data() + offset, read_size);
    offset += read_size;
}


SGEGBox::SGEGBox(const SGEGBox& other)
    : revision(other.revision)
    , spectral_groups(other.spectral_groups)
    , gray_groups(other.gray_groups)
    , exr_attributes(other.exr_attributes)
{}


SGEGBox::~SGEGBox() {}


void SGEGBox::getRaw(std::vector<uint8_t> &data) const
{
    data.resize(size());

    getRaw(data.data());
}

size_t SGEGBox::getRaw(uint8_t data[]) const
{
    size_t offset = 0, write_size = 0;
    uint32_t data_len = 0;

    // revision
    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &revision, write_size);
    offset += write_size;
    
    // spectral_groups
    data_len = spectral_groups.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    for (const SGEGSpectralGroup& sg: spectral_groups) {
        offset += sg.getRaw(data + offset);
    }

    // gray_groups
    data_len = gray_groups.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    for (const SGEGGrayGroup& gg: gray_groups) {
        offset += gg.getRaw(data + offset);
    }

    // exr attributes
    data_len = exr_attributes.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(uint8_t);
    std::memcpy(data + offset, exr_attributes.data(), write_size);
    offset += write_size;

    return offset;
}


size_t SGEGBox::size() const {
    size_t sz = 0;

    sz += sizeof(uint32_t) // + main_layer_name.size() * sizeof(char);
        + 2 * sizeof(uint32_t); // Arrays size;
        
    for (const SGEGSpectralGroup& sg: spectral_groups) {
        sz += sg.size();
    }
    
    for (const SGEGGrayGroup& gg: gray_groups) {
        sz += gg.size();
    }

    // exr_attributes
    sz += sizeof(uint32_t) // Array size
        + exr_attributes.size() * sizeof(uint8_t);

    return sz;
}


SGEGBox& SGEGBox::operator=(const SGEGBox& other)
{        
    if (this != &other) {
        revision        = other.revision;
        spectral_groups = other.spectral_groups;
        gray_groups     = other.gray_groups;
        exr_attributes  = other.exr_attributes;
    }

    return *this;
}

#include <iostream>

void SGEGBox::print() const {
    std::cout << "revision: " << revision << std::endl;
    // std::cout << "main_layer_name: " << main_layer_name.data() << std::endl;
    std::cout << "spectral_groups: " << spectral_groups.size() << std::endl;

    for (const SGEGSpectralGroup& sg: spectral_groups) {
        std::cout << "    root_name: " << sg.root_name.data() << std::endl;

        std::cout << "    layer_indices: ";

        for (const auto& i: sg.layer_indices) {
            std::cout << i << ", ";
        }

        std::cout << std::endl;

        std::cout << "    wavelengths: ";

        for (const auto& i: sg.wavelengths) {
            std::cout << i << ", ";
        }

        std::cout << std::endl;

        std::cout << "    mins: ";

        for (const auto& i: sg.mins) {
            std::cout << i << ", ";
        }

        std::cout << std::endl;

        std::cout << "    maxs: ";

        for (const auto& i: sg.maxs) {
            std::cout << i << ", ";
        }

        std::cout << std::endl << std::endl;        
    }

    std::cout << "gray_groups: " << gray_groups.size() << std::endl;

    for (const SGEGGrayGroup& gg: gray_groups) {
        std::cout << "    layer_name:  "<< gg.layer_name.data() << std::endl;
        std::cout << "    layer_index: "<< gg.layer_index << std::endl;
    }
}
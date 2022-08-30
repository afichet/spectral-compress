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


size_t SGEGGrayGroup::size() const
{
    return
        sizeof(uint32_t) // Array size
        + layer_name.size() * sizeof(char)
        + sizeof(uint32_t);
}


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
        // root_name
        read_size = sizeof(uint32_t);
        std::memcpy(&data_len, data.data() + offset, read_size);
        offset += read_size;

        spectral_groups[i].root_name.resize(data_len);

        read_size = data_len * sizeof(char);
        std::memcpy(spectral_groups[i].root_name.data(), data.data() + offset, read_size);
        offset += read_size;

        // layer_indices
        read_size = sizeof(uint32_t);
        std::memcpy(&data_len, data.data() + offset, read_size);
        offset += read_size;

        spectral_groups[i].layer_indices.resize(data_len);

        read_size = data_len * sizeof(uint32_t);
        std::memcpy(spectral_groups[i].layer_indices.data(), data.data() + offset, read_size);
        offset += read_size;

        // wavelengths
        read_size = sizeof(uint32_t);
        std::memcpy(&data_len, data.data() + offset, read_size);
        offset += read_size;

        spectral_groups[i].wavelengths.resize(data_len);

        read_size = data_len * sizeof(float);
        std::memcpy(spectral_groups[i].wavelengths.data(), data.data() + offset, read_size);
        offset += read_size;

        // mins / maxs
        read_size = sizeof(uint32_t);
        std::memcpy(&data_len, data.data() + offset, read_size);
        offset += read_size;

        spectral_groups[i].mins.resize(data_len);
        spectral_groups[i].maxs.resize(data_len);

        read_size = data_len * sizeof(float);
        std::memcpy(spectral_groups[i].mins.data(), data.data() + offset, read_size);
        offset += read_size;
        std::memcpy(spectral_groups[i].maxs.data(), data.data() + offset, read_size);
        offset += read_size;
    } 

    // gray_groups
    read_size = sizeof(uint32_t);
    std::memcpy(&data_len, data.data() + offset, read_size);
    offset += read_size;

    gray_groups.resize(data_len);

    for (size_t i = 0; i < gray_groups.size(); i++) {
        // layer_name
        read_size = sizeof(uint32_t);
        std::memcpy(&data_len, data.data() + offset, read_size);
        offset += read_size;

        gray_groups[i].layer_name.resize(data_len);

        read_size = data_len * sizeof(char);
        std::memcpy(gray_groups[i].layer_name.data(), data.data() + offset, read_size);
        offset += read_size;

        // layer_index
        read_size = sizeof(uint32_t);
        std::memcpy(&gray_groups[i].layer_index, data.data() + offset, read_size);
        offset += read_size;
    }
}


SGEGBox::SGEGBox(const SGEGBox& other)
    : revision(other.revision)
    , spectral_groups(other.spectral_groups.size())
    , gray_groups(other.gray_groups.size())
{
    for (size_t i = 0; i < spectral_groups.size(); i++) {
        spectral_groups[i].root_name     = other.spectral_groups[i].root_name;
        spectral_groups[i].layer_indices = other.spectral_groups[i].layer_indices;
        spectral_groups[i].wavelengths   = other.spectral_groups[i].wavelengths;
        spectral_groups[i].mins          = other.spectral_groups[i].mins;
        spectral_groups[i].maxs          = other.spectral_groups[i].maxs;
    }

    for (size_t i = 0; i < gray_groups.size(); i++) {
        gray_groups[i].layer_name  = other.gray_groups[i].layer_name;
        gray_groups[i].layer_index = other.gray_groups[i].layer_index;
    }
}


SGEGBox::~SGEGBox() {}


void SGEGBox::getRaw(std::vector<uint8_t> &data) const
{
    data.resize(size());

    size_t offset = 0, write_size = 0;
    uint32_t data_len = 0;

    // revision
    write_size = sizeof(uint32_t);
    std::memcpy(data.data() + offset, &revision, write_size);
    offset += write_size;
    
    // spectral_groups
    data_len = spectral_groups.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data.data() + offset, &data_len, write_size);
    offset += write_size;

    for (const SGEGSpectralGroup& sg: spectral_groups) {
        // root_name
        data_len = sg.root_name.size();

        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &data_len, write_size);
        offset += write_size;

        write_size = data_len * sizeof(char);
        std::memcpy(data.data() + offset, sg.root_name.data(), write_size);
        offset += write_size;

        // layer_indices
        data_len = sg.layer_indices.size();

        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &data_len, write_size);
        offset += write_size;

        write_size = data_len * sizeof(uint32_t);
        std::memcpy(data.data() + offset, sg.layer_indices.data(), write_size);
        offset += write_size;

        // wavelengths
        data_len = sg.wavelengths.size();
        
        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &data_len, write_size);
        offset += write_size;

        write_size = data_len * sizeof(float);
        std::memcpy(data.data() + offset, sg.wavelengths.data(), write_size);
        offset += write_size;

        // mins / maxs
        data_len = sg.mins.size();

        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &data_len, write_size);
        offset += write_size;

        write_size = data_len * sizeof(float);
        std::memcpy(data.data() + offset, sg.mins.data(), write_size);
        offset += write_size;
        std::memcpy(data.data() + offset, sg.maxs.data(), write_size);
        offset += write_size;
    }

    // gray_groups
    data_len = gray_groups.size();

    write_size = sizeof(uint32_t);
    std::memcpy(data.data() + offset, &data_len, write_size);
    offset += write_size;

    for (const SGEGGrayGroup& gg: gray_groups) {
        // layer_name
        data_len = gg.layer_name.size();

        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &data_len, write_size);
        offset += write_size;

        write_size = data_len * sizeof(char);
        std::memcpy(data.data() + offset, gg.layer_name.data(), write_size);
        offset += write_size;

        // layer_index
        write_size = sizeof(uint32_t);
        std::memcpy(data.data() + offset, &gg.layer_index, write_size);
        offset += write_size;
    }
}


size_t SGEGBox::size() const {
    size_t sz = sizeof(uint32_t);// + main_layer_name.size() * sizeof(char);
    sz += 2 * sizeof(uint32_t); // Arrays size;
        
    for (const SGEGSpectralGroup& sg: spectral_groups) {
        sz += sg.size();
    }
    
    for (const SGEGGrayGroup& gg: gray_groups) {
        sz += gg.size();
    }

    return sz;
}


SGEGBox& SGEGBox::operator=(const SGEGBox& other)
{        
    if (this != &other) {
        revision = other.revision;
        // main_layer_name = other.main_layer_name;

        spectral_groups.resize(other.spectral_groups.size());

        for (size_t i = 0; i < spectral_groups.size(); i++) {
            spectral_groups[i].root_name     = other.spectral_groups[i].root_name;
            spectral_groups[i].layer_indices = other.spectral_groups[i].layer_indices;
            spectral_groups[i].wavelengths   = other.spectral_groups[i].wavelengths;
            spectral_groups[i].mins          = other.spectral_groups[i].mins;
            spectral_groups[i].maxs          = other.spectral_groups[i].maxs;
        }

        gray_groups.resize(other.gray_groups.size());

        for (size_t i = 0; i < gray_groups.size(); i++) {
            gray_groups[i].layer_name  = other.gray_groups[i].layer_name;
            gray_groups[i].layer_index = other.gray_groups[i].layer_index;
        }
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
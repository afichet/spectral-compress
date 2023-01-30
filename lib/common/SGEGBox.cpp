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

#include "SGEGBox.h"

#include <cstring>
#include <cstdint>

template<typename T>
size_t write_value(const T& value, uint8_t data[])
{
    const size_t write_size = sizeof(T);

    std::memcpy(data, &value, write_size);

    return write_size;
}

template<typename T>
size_t read_value(const uint8_t data[], T& value)
{
    const size_t read_size = sizeof(T);

    std::memcpy(&value, data, read_size);

    return read_size;
}

template<typename T>
size_t write_vector(const std::vector<T>& vec, uint8_t data[])
{
    size_t offset = 0, write_size = 0;
    uint32_t data_len = 0;

    data_len = vec.size();

    write_size = sizeof(data_len);
    std::memcpy(data + offset, &data_len, write_size);
    offset += write_size;

    write_size = data_len * sizeof(T);
    std::memcpy(data + offset, vec.data(), write_size);
    offset += write_size;

    return offset;
}

template<typename T>
size_t read_vector(const uint8_t data[], std::vector<T>& vec)
{
    size_t offset = 0, read_size = 0;
    uint32_t data_len = 0;

    read_size = sizeof(data_len);
    std::memcpy(&data_len, data + offset, read_size);
    offset += read_size;

    vec.resize(data_len);

    read_size = data_len * sizeof(T);
    std::memcpy(vec.data(), data + offset, read_size);
    offset += read_size;

    return offset;
}

// ----------------------------------------------------------------------------

SGEGSpectralGroup::SGEGSpectralGroup() {}


SGEGSpectralGroup::SGEGSpectralGroup(const SGEGSpectralGroup& other)
    : root_name(other.root_name)
    , layer_indices(other.layer_indices)
    , wavelengths(other.wavelengths)
    , mins(other.mins)
    , maxs(other.maxs)
    , method(other.method)
{
}


size_t SGEGSpectralGroup::getRaw(uint8_t data[]) const
{
    size_t offset = 0;

    offset += write_vector(root_name, data + offset);
    offset += write_vector(layer_indices, data + offset);
    offset += write_vector(wavelengths, data + offset);
    offset += write_vector(mins, data + offset);
    offset += write_vector(maxs, data + offset);
    offset += write_value(method, data + offset);

    return offset;
}


size_t SGEGSpectralGroup::fromRaw(const uint8_t data[])
{
    size_t offset = 0;

    offset += read_vector(data + offset, root_name);
    offset += read_vector(data + offset, layer_indices);
    offset += read_vector(data + offset, wavelengths);
    offset += read_vector(data + offset, mins);
    offset += read_vector(data + offset, maxs);
    offset += read_value(data + offset, method);

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
        + sizeof(uint32_t) // Arrays size
        + maxs.size() * sizeof(float)
        + sizeof(float)
        + sizeof(float)
        + sizeof(SpectralCompressionType);
}

// ----------------------------------------------------------------------------

SGEGGrayGroup::SGEGGrayGroup() {}


SGEGGrayGroup::SGEGGrayGroup(const SGEGGrayGroup& other)
    : layer_name(other.layer_name)
    , layer_index(other.layer_index)
{
}

size_t SGEGGrayGroup::getRaw(uint8_t data[]) const
{
    size_t offset = 0;

    offset += write_vector(layer_name, data + offset);
    offset += write_value(layer_index, data + offset);

    return offset;
}


size_t SGEGGrayGroup::fromRaw(const uint8_t data[])
{
    size_t offset = 0;

    offset += read_vector(data + offset, layer_name);
    offset += read_value(data + offset, layer_index);

    return offset;
}


size_t SGEGGrayGroup::size() const
{
    return
        sizeof(uint32_t) // Array size
        + layer_name.size() * sizeof(char)
        + sizeof(uint32_t);
}

// ----------------------------------------------------------------------------

SGEGBox::SGEGBox()
    : revision(1)
    , n_parts(1)
{}


SGEGBox::SGEGBox(const std::vector<uint8_t> &data)
{
    size_t offset = 0;
    uint32_t n_spectral_groups, n_gray_groups;

    offset += read_value(data.data() + offset, revision);
    offset += read_value(data.data() + offset, n_parts);

    // spectral groups
    offset += read_value(data.data() + offset, n_spectral_groups);

    spectral_groups.resize(n_spectral_groups);

    for (size_t i = 0; i < spectral_groups.size(); i++) {
        offset += spectral_groups[i].fromRaw(data.data() + offset);
    }

    // gray_groups
    offset += read_value(data.data() + offset, n_gray_groups);

    gray_groups.resize(n_gray_groups);

    for (size_t i = 0; i < gray_groups.size(); i++) {
        offset += gray_groups[i].fromRaw(data.data() + offset);
    }

    // exr attributes
    offset += read_vector(data.data() + offset, exr_attributes);
}


SGEGBox::SGEGBox(const SGEGBox& other)
    : revision(other.revision)
    , n_parts(other.n_parts)
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
    size_t offset = 0;
    const uint32_t n_spectral_groups = spectral_groups.size();
    const uint32_t n_gray_groups = gray_groups.size();

    // revision
    offset += write_value(revision, data + offset);

    // number of parts
    offset += write_value(n_parts, data + offset);

    // spectral_groups
    offset += write_value(n_spectral_groups, data + offset);

    for (const SGEGSpectralGroup& sg: spectral_groups) {
        offset += sg.getRaw(data + offset);
    }

    // gray_groups
    offset += write_value(n_gray_groups, data + offset);

    for (const SGEGGrayGroup& gg: gray_groups) {
        offset += gg.getRaw(data + offset);
    }

    // exr attributes
    offset += write_vector(exr_attributes, data + offset);

    return offset;
}


size_t SGEGBox::size() const {
    size_t sz = 0;

    sz += sizeof(revision);
    sz += sizeof(n_parts);

    sz += sizeof(uint32_t); // spectral_group array size

    for (const SGEGSpectralGroup& sg: spectral_groups) {
        sz += sg.size();
    }

    sz += sizeof(uint32_t); // gray_groups array size

    for (const SGEGGrayGroup& gg: gray_groups) {
        sz += gg.size();
    }

    // exr_attributes
    sz += sizeof(uint32_t); // Array size
    sz += exr_attributes.size() * sizeof(uint8_t);

    return sz;
}


SGEGBox& SGEGBox::operator=(const SGEGBox& other)
{
    if (this != &other) {
        revision        = other.revision;
        n_parts         = other.n_parts;
        spectral_groups = other.spectral_groups;
        gray_groups     = other.gray_groups;
        exr_attributes  = other.exr_attributes;
    }

    return *this;
}

#include <iostream>

void SGEGBox::print() const {
    std::cout << "revision: " << revision << std::endl;
    std::cout << "n_parts: " << n_parts << std::endl;
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

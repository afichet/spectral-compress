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

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <string>

#include <spectral_compression_type.h>


struct SGEGSpectralGroup
{
    SGEGSpectralGroup();
    SGEGSpectralGroup(const SGEGSpectralGroup& other);

    // Contains everything before the last `.` separating the layer name from
    // the spectral wavelength.
    // E.g., `T.480`       has `T` as root name,
    //       `left.S0.340` has `left.S0` as root name.
    std::vector<char> root_name;
    // Corresponding indices of the framebuffer in the JXL file
    // If the file is split accros multiple files (when the number of
    // framebuffers exceed the 256 limit), this value can exceed 256 and is
    // refering to the n_th layer as if we've used a single file.
    std::vector<uint32_t> layer_indices;

    // The original wavelengths the signal was sampled at
    std::vector<float> wavelengths;

    // Mins & Maxs of the moments starting from 1st moment,
    // 0th moment is not rescaled
    std::vector<float> mins;
    std::vector<float> maxs;

    SpectralCompressionType method;

    /**
     * Get the content of this structure for serialization.
     * @param data An array of bytes which need to be at least of `size`
     *             elements.
     * @return The number of bytes written
     */
    size_t getRaw(uint8_t data[]) const;

    size_t fromRaw(const uint8_t data[]);

    /**
     * Size in bytes require to serialize this data structure i.e., the minimum
     * array size required when calling `getRaw`
     */
    size_t size() const;

    // SGEGSpectralGroup& operator=(const SGEGSpectralGroup& other);
};


struct SGEGGrayGroup
{
    SGEGGrayGroup();
    SGEGGrayGroup(const SGEGGrayGroup& other);

    // The name the layer has in the original OpenEXR file
    std::vector<char> layer_name;
    // Corresponding index of the framebuffer in the JXL file.
    // If the file is split accros multiple files (when the number of
    // framebuffers exceed the 256 limit), this value can exceed 256 and is
    // refering to the n_th layer as if we've used a single file.
    uint32_t layer_index;

    /**
     * Get the content of this structure for serialization.
     * @param data An array of bytes which need to be at least of `size`
     *             elements.
     * @return The number of bytes written
     */
    size_t getRaw(uint8_t data[]) const;

    size_t fromRaw(const uint8_t data[]);

    /**
     * Size in bytes require to serialize this data structure i.e., the minimum
     * array size required when calling `getRaw`
     */
    size_t size() const;
};


// Spectral Graphics Extended Group

struct SGEGBox
{
    // Revision of the format, curretnly at 1.0
    uint32_t revision;
    // Number of files needed to store this image: at the time of writing this
    // code, JXL does not support more than 1 main image and 255 subimages
    uint32_t n_parts;

    // List of images, either gray ones (a single framebuffer) or spectral ones
    std::vector<SGEGSpectralGroup> spectral_groups;
    std::vector<SGEGGrayGroup> gray_groups;

    // Raw OpenEXR attributes from the original OpenEXR image.
    std::vector<uint8_t> exr_attributes;

    SGEGBox();

    SGEGBox(const std::vector<uint8_t> &data);

    SGEGBox(const SGEGBox& other);

    virtual ~SGEGBox();

    /**
     * Get the content of this structure for serialization.
     * @param data An array of bytes
     * @return The number of bytes written
     */
    void getRaw(std::vector<uint8_t> &data) const;

    /**
     * Get the content of this structure for serialization.
     * @param data An array of bytes which need to be at least of `size`
     *             elements.
     * @return The number of bytes written
     */
    size_t getRaw(uint8_t data[]) const;

    /**
     * Size in bytes require to serialize this data structure i.e., the minimum
     * array size required when calling `getRaw`
     */
    size_t size() const;

    SGEGBox& operator=(const SGEGBox& other);

    /**
     * Pretty print on `stdout` the content of this structure
     */
    void print() const;
};

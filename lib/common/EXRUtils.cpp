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

#include "EXRUtils.h"


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


void vec_float_to_uint32(
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


char* condition_framebuffers(
    const std::vector<float>& framebuffer_in,
    Imf::PixelType write_pixel_type,
    size_t& type_stride,
    std::vector<std::vector<uint32_t>>& framebuffers_int,
    std::vector<std::vector<half>>& framebuffers_half)
{
    char* fb_data;

    // TODO: Natively support uint32 without a loss of precision when the data
    // are provided in this form
    switch (write_pixel_type)
    {
        case Imf::PixelType::UINT:
            type_stride = sizeof(uint32_t);
            framebuffers_int.resize(framebuffers_int.size() + 1);
            vec_float_to_uint32(framebuffer_in, framebuffers_int.back());
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

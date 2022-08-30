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

#include "JXLImage.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>

// ----------------------------------------------------------------------------

#define CHECK_JXL_ENC_STATUS(status)                                          \
    if (JXL_ENC_SUCCESS != status) {                                          \
        std::cerr << "[ERR] At: "                                             \
                  << __FILE__  << ":" << __LINE__                             \
                  << " " << status << std::endl;                              \
    }                                                                         \


#define CHECK_JXL_DEC_STATUS(status)                                          \
    if (JXL_DEC_SUCCESS != status) {                                          \
        std::cerr << "[ERR] At: "                                             \
                  << __FILE__  << ":" << __LINE__                             \
                  << " " << status << std::endl;                              \
    }                                                                         \

// ----------------------------------------------------------------------------

JXLFramebuffer::JXLFramebuffer(
    uint32_t width, uint32_t height,
    uint32_t n_color_channels,
    uint32_t n_bits_per_sample, 
    uint32_t n_exponent_bits_per_sample,
    uint32_t downsampling_factor,
    const char* name)
    : _name(nullptr)
    , _n_bits_per_sample(n_bits_per_sample)
    , _n_exponent_bits_per_sample(n_exponent_bits_per_sample)
    , _downsampling_factor(downsampling_factor)
    , _pixel_format(
        {
            /* .num_channels = */ n_color_channels,
            /* .data_type    = */ JXL_TYPE_FLOAT,
            /* .endianness   = */ JXL_NATIVE_ENDIAN,
            /* .align        = */ 0
        })
    , _pixel_data(width * height * n_color_channels)
{
    if (name != nullptr) { setName(name); }
}



JXLFramebuffer::JXLFramebuffer(
    const std::vector<float>& framebuffer,
    uint32_t n_color_channels,
    uint32_t n_bits_per_sample, 
    uint32_t n_exponent_bits_per_sample,
    uint32_t downsampling_factor,
    const char* name)
    : _name(nullptr)
    , _n_bits_per_sample(n_bits_per_sample)
    , _n_exponent_bits_per_sample(n_exponent_bits_per_sample)
    , _downsampling_factor(downsampling_factor)
    , _pixel_format(
        {
            /* .num_channels = */ n_color_channels,
            /* .data_type    = */ JXL_TYPE_FLOAT,
            /* .endianness   = */ JXL_NATIVE_ENDIAN,
            /* .align        = */ 0
        })
    , _pixel_data(framebuffer)
{
    if (name != nullptr) { setName(name); }
}


JXLFramebuffer::~JXLFramebuffer()
{
    delete[] _name;
}


void JXLFramebuffer::setName(const char* name) {
    delete[] _name;
    _name = nullptr;

    if (name != nullptr) {
        const size_t len_name = std::strlen(name);
        _name = new char[len_name + 1];
        std::strcpy(_name, name);
    }
}


// ===========================================================================

// TODO: remove
#include <iostream>

JXLImage::JXLImage(uint32_t width, uint32_t height)
    : _width(width)
    , _height(height)
{}


JXLImage::JXLImage(const char* filename)
{
    JxlDecoderStatus status;
    
    std::vector<uint8_t> jxl_data;

    JxlDecoderPtr              dec;
    JxlThreadParallelRunnerPtr runner;
    std::vector<uint8_t> box_raw_sgeg;

    // Read file
    FILE* file = fopen(filename, "rb");

    fseek(file, 0, SEEK_END);
    const long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    jxl_data.resize(file_size);

    fread(jxl_data.data(), 1, file_size, file);
    fclose(file);

    // Decode JXL stream
    dec = JxlDecoderMake(nullptr);
    runner = JxlThreadParallelRunnerMake(nullptr, JxlThreadParallelRunnerDefaultNumWorkerThreads());

    status = JxlDecoderSetParallelRunner(
        dec.get(), 
        JxlThreadParallelRunner, 
        runner.get()
    );

    CHECK_JXL_DEC_STATUS(status);

    status = JxlDecoderSubscribeEvents(
        dec.get(), 
        JXL_DEC_BASIC_INFO |
        JXL_DEC_BOX |
        // JXL_DEC_COLOR_ENCODING | 
        JXL_DEC_FULL_IMAGE);
    
    CHECK_JXL_DEC_STATUS(status);

    JxlDecoderSetInput(dec.get(), jxl_data.data(), jxl_data.size());

    for (JxlDecoderStatus status_process = JXL_DEC_NEED_MORE_INPUT; 
        status_process != JXL_DEC_FULL_IMAGE; 
        status_process = JxlDecoderProcessInput(dec.get())) 
    {
        switch (status_process) {
            case JXL_DEC_BASIC_INFO:
                {
                    // Read metadata and allocate memory
                    JxlBasicInfo basic_info;
                    JxlExtraChannelInfo extra_info;
                    std::vector<char> layer_name;

                    // Main layer metadata
                    status = JxlDecoderGetBasicInfo(dec.get(), &basic_info);

                    _width  = basic_info.xsize;
                    _height = basic_info.ysize;

                    _framebuffers.push_back(new JXLFramebuffer(
                        _width, _height,
                        basic_info.num_color_channels,
                        basic_info.bits_per_sample,
                        basic_info.exponent_bits_per_sample
                    ));

                    // Extra layer metadata
                    for (size_t i = 0; i < basic_info.num_extra_channels; i++) {
                        status = JxlDecoderGetExtraChannelInfo(dec.get(), i, &extra_info);
                        CHECK_JXL_DEC_STATUS(status);

                        if (extra_info.name_length > 0) {
                            layer_name.resize(extra_info.name_length + 1);
                            status = JxlDecoderGetExtraChannelName(dec.get(), i, layer_name.data(), extra_info.name_length);
                            layer_name[extra_info.name_length] = 0;
                        }

                        _framebuffers.push_back(new JXLFramebuffer(
                            _width, _height,
                            basic_info.num_color_channels, // TODO, is this really the case?
                            extra_info.bits_per_sample,
                            extra_info.exponent_bits_per_sample,
                            1,
                            (extra_info.name_length > 0) ? layer_name.data() : nullptr
                        ));
                    }
                }
                break;

            case JXL_DEC_NEED_IMAGE_OUT_BUFFER:
                {
                    JxlPixelFormat format = _framebuffers[0]->getPixelFormat();

                    status = JxlDecoderSetImageOutBuffer(
                        dec.get(), 
                        &format,
                        _framebuffers[0]->getPixelData().data(),
                        _framebuffers[0]->getSizeBytes()
                    );
                    CHECK_JXL_DEC_STATUS(status);

                    for (size_t i = 1; i < _framebuffers.size(); i++) {
                        format = _framebuffers[i]->getPixelFormat();

                        status = JxlDecoderSetExtraChannelBuffer(
                            dec.get(),
                            &format,
                            _framebuffers[i]->getPixelData().data(),
                            _framebuffers[i]->getSizeBytes(),
                            i - 1
                        );
                        CHECK_JXL_DEC_STATUS(status);
                    }
                }
                break;

            case JXL_DEC_BOX:
                {
                    JxlBoxType box_type;
                    
                    status = JxlDecoderGetBoxType(dec.get(), box_type, JXL_TRUE);
                    CHECK_JXL_DEC_STATUS(status);

                    if (box_type[0] == 's' 
                     && box_type[1] == 'g' 
                     && box_type[2] == 'e' 
                     && box_type[3] == 'g') {
                        uint64_t box_size;
                        status = JxlDecoderGetBoxSizeRaw(dec.get(), &box_size);
                        box_raw_sgeg.resize(box_size);

                        JxlDecoderSetBoxBuffer(dec.get(), box_raw_sgeg.data(), box_raw_sgeg.size());
                    }
                }
                break;
            case JXL_DEC_NEED_MORE_INPUT:
                break;
            case JXL_DEC_ERROR:
                CHECK_JXL_DEC_STATUS(status_process);
                break;
            default:
                std::cout << "Unknown decoder status: "
                          << status_process << std::endl;
                break;
        }
    }

    JxlDecoderReleaseInput(dec.get());
    
    _sgeg_box = SGEGBox(box_raw_sgeg);

    CHECK_JXL_DEC_STATUS(status);
}

JXLImage::~JXLImage()
{
    for (size_t i = 0; i < _framebuffers.size(); i++) {
        delete _framebuffers[i];
    }
}


void JXLImage::setBox(const SGEGBox& box)
{
    _sgeg_box = box;
}


size_t JXLImage::appendFramebuffer(
    const std::vector<float>& framebuffer,
    uint32_t n_channels,
    uint32_t enc_bits_per_sample,
    uint32_t enc_exponent_bits_per_sample,
    uint32_t enc_downsampling_factor,
    const char* name)
{
    assert(framebuffer.size() == n_channels * _width * _height);

    JXLFramebuffer *fb = new JXLFramebuffer(
        framebuffer, 
        n_channels,
        enc_bits_per_sample,
        enc_exponent_bits_per_sample,
        enc_downsampling_factor,
        name
    );

    _framebuffers.push_back(fb);

    return _framebuffers.size() - 1;
}


void JXLImage::write(const char* filename) const {
    // TODO: Check if _framebuffers.size() > 0
    JxlEncoderStatus           status;
    JxlEncoderPtr              enc;
    JxlThreadParallelRunnerPtr runner;
    JxlBasicInfo               basic_info;
    JxlColorEncoding           color_encoding;
    
    runner = JxlThreadParallelRunnerMake(
            nullptr, 
            JxlThreadParallelRunnerDefaultNumWorkerThreads()
        );

    enc = JxlEncoderMake(nullptr);

    status = JxlEncoderSetParallelRunner(
        enc.get(),
        JxlThreadParallelRunner,
        runner.get()
    );

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------

    status = JxlEncoderUseBoxes(enc.get());

    CHECK_JXL_ENC_STATUS(status);

    char tp[4] = {'s','g','e','g'};

    std::vector<uint8_t> raw_box;
    _sgeg_box.getRaw(raw_box);
    
    status = JxlEncoderAddBox(enc.get(), tp, raw_box.data(), raw_box.size(), JXL_FALSE);
    
    CHECK_JXL_ENC_STATUS(status);

    JxlEncoderCloseBoxes(enc.get());

    // ====================================================================
    // Main framebuffer
    // ====================================================================

    JxlEncoderInitBasicInfo(&basic_info);

    basic_info.xsize                    = _width;
    basic_info.ysize                    = _height;
    basic_info.num_extra_channels       = _framebuffers.size() - 1;
    basic_info.num_color_channels       = _framebuffers[0]->getNColorChannels();
    basic_info.bits_per_sample          = _framebuffers[0]->getBitsPerSample();
    basic_info.exponent_bits_per_sample = _framebuffers[0]->getExponentBitsPerSample();
    basic_info.uses_original_profile    = JXL_TRUE;

    status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------

    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

    // Set compression quality
    // status = JxlEncoderSetFrameDistance(frame_settings, .1);
    // CHECK_JXL_ENC_STATUS(status);

    // status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    // CHECK_JXL_DEC_STATUS(status);

    // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
    status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, _framebuffers[0]->getDownsamplingFactor());
    CHECK_JXL_ENC_STATUS(status);

    // status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    // CHECK_JXL_ENC_STATUS(status);

    status = JxlEncoderSetFrameDistance(frame_settings, 0.);
    CHECK_JXL_ENC_STATUS(status);

    // TODO: Fixme!
    switch (_framebuffers[0]->getNColorChannels()) {
        case 1:
            color_encoding.color_space = JXL_COLOR_SPACE_GRAY;
            break;
        case 3:
        case 4:
            color_encoding.color_space = JXL_COLOR_SPACE_RGB;
            break;
        default:
            std::cerr << "Unknown color space matching number of channels" << std::endl;
            assert(0);
            break;
    }

    color_encoding.white_point       = JXL_WHITE_POINT_D65;
    color_encoding.primaries         = JXL_PRIMARIES_SRGB;
    color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
    color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

    // TODO set name

    status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);

    CHECK_JXL_ENC_STATUS(status);

    const JxlPixelFormat format = _framebuffers[0]->getPixelFormat();

    const size_t data_size = _width * _height * _framebuffers[0]->getNColorChannels() * sizeof(float);

    status = JxlEncoderAddImageFrame(
        frame_settings, 
        &format,
        _framebuffers[0]->getPixelData().data(),
        data_size
    );

    CHECK_JXL_ENC_STATUS(status);

    // ====================================================================
    // Subframebuffers
    // ====================================================================

    for (size_t i = 1; i < _framebuffers.size(); i++) {
        const JXLFramebuffer* fb = _framebuffers[i];
        const JxlPixelFormat format = fb->getPixelFormat();

        JxlExtraChannelInfo extra_info;
        JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);

        extra_info.bits_per_sample          = fb->getBitsPerSample();
        extra_info.exponent_bits_per_sample = fb->getExponentBitsPerSample();

        JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

        // Set compression quality
        // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
        // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EXTRA_CHANNEL_RESAMPLING, downsampling);
        // status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EXTRA_CHANNEL_RESAMPLING, 4);
        // CHECK_JXL_ENC_STATUS(status);

        // // TODO: Fixme!
        // color_encoding.color_space       = JXL_COLOR_SPACE_GRAY;
        // color_encoding.white_point       = JXL_WHITE_POINT_D65;
        // color_encoding.primaries         = JXL_PRIMARIES_SRGB;
        // color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
        // color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

        // status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);

        const size_t data_size = _width * _height * fb->getNColorChannels() * sizeof(float);

        status = JxlEncoderSetExtraChannelInfo(enc.get(), i - 1, &extra_info);
        CHECK_JXL_ENC_STATUS(status);

        if (fb->getName() != nullptr) {
            const char* channel_name = fb->getName();
            status = JxlEncoderSetExtraChannelName(enc.get(), 0, channel_name, std::strlen(channel_name));
            CHECK_JXL_ENC_STATUS(status);
        }

        status = JxlEncoderSetExtraChannelBuffer(
            frame_settings, 
            &format, 
            fb->getPixelDataConst().data(),
            data_size,
            i - 1);

        CHECK_JXL_ENC_STATUS(status);
    }

    // ------------------------------------------------------------------------

    JxlEncoderCloseInput(enc.get());

    std::vector<uint8_t> compressed(64);

    uint8_t* next_out = compressed.data();
    size_t avail_out = compressed.size() - (next_out - compressed.data());
    
    status = JXL_ENC_NEED_MORE_OUTPUT;

    // for (status = JXL_ENC_NEED_MORE_OUTPUT ;; status == JXL_ENC_NEED_MORE_OUTPUT) {
    while (status == JXL_ENC_NEED_MORE_OUTPUT) {
        status = JxlEncoderProcessOutput(enc.get(), &next_out, &avail_out);

        if (status == JXL_ENC_NEED_MORE_OUTPUT) {
            size_t offset = next_out - compressed.data();
            compressed.resize(compressed.size() * 2);

            next_out = compressed.data() + offset;
            avail_out = compressed.size() - offset;
        }
    }

    CHECK_JXL_ENC_STATUS(status);
    
    compressed.resize(next_out - compressed.data());

    // Write file
    FILE* file = fopen(filename, "wb");
    fwrite(compressed.data(), sizeof(uint8_t), compressed.size(), file);
    fclose(file);
}


// void JXLImageReader::print_basic_info() const
// {
//     std::cout << "┌                   width: " << _basic_info.xsize << std::endl;
//     std::cout << "├                  height: " << _basic_info.ysize << std::endl;
//     std::cout << "├         bits per sample: " << _basic_info.bits_per_sample << std::endl;
//     std::cout << "├exponent bits per sample: " << _basic_info.exponent_bits_per_sample << std::endl;
//     std::cout << "├        intensity target: " << _basic_info.intensity_target << std::endl;
//     std::cout << "├ relative to max display: " << _basic_info.relative_to_max_display << std::endl;
//     std::cout << "├            linear below: " << _basic_info.linear_below << std::endl;
//     std::cout << "├   uses original profile: " << _basic_info.uses_original_profile << std::endl;
//     std::cout << "├            have preview: " << _basic_info.have_preview << std::endl;
//     std::cout << "├          have animation: " << _basic_info.have_animation << std::endl;
//     std::cout << "├             orientation: " << _basic_info.orientation << std::endl;
//     std::cout << "├      num color channels: " << _basic_info.num_color_channels << std::endl;
//     std::cout << "├      num extra channels: " << _basic_info.num_extra_channels << std::endl;
//     std::cout << "├              alpha bits: " << _basic_info.alpha_bits << std::endl;
//     std::cout << "├     alpha exponent bits: " << _basic_info.alpha_exponent_bits << std::endl;
//     std::cout << "└     alpha premultiplied: " << _basic_info.alpha_premultiplied << std::endl;


//     for (uint32_t i = 0; i < _basic_info.num_extra_channels; i++) {
//         std::cout << "\t┌ Channel #" << i << std::endl;
//         std::cout << "\t└  name: " << _sub_framebuffers_names[i] << std::endl; 
//     }

//     for (uint32_t i = 0; i < _sgeg_box.n_moments; i++) {
//         std::cout << "xmin: " << _sgeg_box.moment_min[i] << std::endl;
//         std::cout << "xmax: " << _sgeg_box.moment_max[i] << std::endl;
//     }
// }
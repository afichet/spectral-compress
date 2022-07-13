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


JXLImageWriter::JXLImageWriter(
    size_t width,
    size_t height,
    SGEG_box sgeg_box,
    size_t n_sub_framebuffers)
    : _enc(JxlEncoderMake(nullptr))
    , _sgeg_box(sgeg_box)
{
    JxlEncoderStatus status;

    // ------------------------------------------------------------------------

    // _enc = JxlEncoderMake(nullptr);
    _runner = JxlThreadParallelRunnerMake(
        nullptr, 
        JxlThreadParallelRunnerDefaultNumWorkerThreads()
    );

    status = JxlEncoderSetParallelRunner(
        _enc.get(),
        JxlThreadParallelRunner,
        _runner.get()
    );

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------

    JxlEncoderInitBasicInfo(&_basic_info);

    _basic_info.xsize                    = width;
    _basic_info.ysize                    = height;
    _basic_info.num_color_channels       = 1;
    _basic_info.num_extra_channels       = n_sub_framebuffers;
    _basic_info.bits_per_sample          = 32;
    _basic_info.exponent_bits_per_sample = 8;
    _basic_info.uses_original_profile    = JXL_TRUE;

    status = JxlEncoderSetBasicInfo(_enc.get(), &_basic_info);

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------

    _color_encoding.color_space       = JXL_COLOR_SPACE_GRAY;
    _color_encoding.white_point       = JXL_WHITE_POINT_D65;
    _color_encoding.primaries         = JXL_PRIMARIES_SRGB;
    _color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
    _color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

    status = JxlEncoderSetColorEncoding(_enc.get(), &_color_encoding);

    CHECK_JXL_ENC_STATUS(status);

    // ------------------------------------------------------------------------
    
    status = JxlEncoderUseBoxes(_enc.get());
    
    CHECK_JXL_ENC_STATUS(status);

    char tp[4] = {'s','g','e','g'};

    std::vector<uint8_t> raw_box;
    _sgeg_box.getRaw(raw_box);
    
    status = JxlEncoderAddBox(_enc.get(), tp, raw_box.data(), raw_box.size(), JXL_FALSE);
    
    CHECK_JXL_ENC_STATUS(status);
}


void JXLImageWriter::addMainFramebuffer(
    void* framebuffer,
    JxlPixelFormat pixel_format,
    size_t el_size,
    int32_t downsampling)
{
    JxlEncoderStatus status;
    
    // Deprecated
    // JxlEncoderOptions* frame_settings = JxlEncoderOptionsCreate(_enc.get(), nullptr);
    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(_enc.get(), nullptr);
    
    // Set compression quality
    // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
    status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, downsampling);
    CHECK_JXL_ENC_STATUS(status);

    status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
    CHECK_JXL_ENC_STATUS(status);

    const size_t data_size = _basic_info.xsize * _basic_info.ysize * el_size;

    status = JxlEncoderAddImageFrame(
        frame_settings, 
        &pixel_format,
        framebuffer,
        data_size
    );

    CHECK_JXL_ENC_STATUS(status);
}


void JXLImageWriter::addMainFramebuffer(float* framebuffer, int32_t downsampling)
{
    JxlPixelFormat format = {1, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};
    addMainFramebuffer(framebuffer, format, sizeof(float), downsampling);
}


void JXLImageWriter::addMainFramebuffer(uint8_t* framebuffer, int32_t downsampling)
{
    JxlPixelFormat format = {1, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};
    addMainFramebuffer(framebuffer, format, sizeof(uint8_t), downsampling);
}


void JXLImageWriter::addMainFramebuffer(uint16_t* framebuffer, int32_t downsampling)
{
    JxlPixelFormat format = {1, JXL_TYPE_UINT16, JXL_NATIVE_ENDIAN, 0};
    addMainFramebuffer(framebuffer, format, sizeof(uint16_t), downsampling);
}


void JXLImageWriter::addSubFramebuffer(
    void* framebuffer,
    JxlPixelFormat pixel_format,
    JxlExtraChannelInfo extra_info,
    size_t el_size,
    size_t index,
    int32_t downsampling,
    const char* channel_name)
{
    JxlEncoderStatus status;

    // Deprecated
    // JxlEncoderOptions* frame_settings = JxlEncoderOptionsCreate(_enc.get(), nullptr);
    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(_enc.get(), nullptr);

    // Set compression quality
    // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, 9);
    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EXTRA_CHANNEL_RESAMPLING, downsampling);

    const size_t data_size = _basic_info.xsize * _basic_info.ysize * el_size;

    status = JxlEncoderSetExtraChannelInfo(_enc.get(), index, &extra_info);
    CHECK_JXL_ENC_STATUS(status);

    if (channel_name) {
        status = JxlEncoderSetExtraChannelName(_enc.get(), 0, channel_name, std::strlen(channel_name));
        CHECK_JXL_ENC_STATUS(status);
    }

    status = JxlEncoderSetExtraChannelBuffer(
        frame_settings, 
        &pixel_format, 
        framebuffer,
        data_size,
        index);

    CHECK_JXL_ENC_STATUS(status);
}


void JXLImageWriter::addSubFramebuffer(
    float* framebuffer,
    size_t index,
    int32_t downsampling,
    const char* channel_name)
{
    JxlPixelFormat format = {1, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};

    JxlExtraChannelInfo extra_info;
    JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);
    extra_info.bits_per_sample = 32;
    extra_info.exponent_bits_per_sample = 8;

    addSubFramebuffer(framebuffer, format, extra_info, sizeof(float), index, downsampling, channel_name);
}


void JXLImageWriter::addSubFramebuffer(
    uint8_t* framebuffer,
    size_t index,
    int32_t downsampling,
    const char* channel_name)
{
    JxlPixelFormat format = {1, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};

    JxlExtraChannelInfo extra_info;
    JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);
    extra_info.bits_per_sample = 8;
    extra_info.exponent_bits_per_sample = 0;

    addSubFramebuffer(framebuffer, format, extra_info, sizeof(uint8_t), index, downsampling, channel_name);
}


void JXLImageWriter::addSubFramebuffer(
    uint16_t* framebuffer,
    size_t index,
    int32_t downsampling,
    const char* channel_name)
{
    JxlPixelFormat format = {1, JXL_TYPE_UINT16, JXL_NATIVE_ENDIAN, 0};

    JxlExtraChannelInfo extra_info;
    JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);
    extra_info.bits_per_sample = 16;
    extra_info.exponent_bits_per_sample = 0;

    addSubFramebuffer(framebuffer, format, extra_info, sizeof(uint16_t), index, downsampling, channel_name);
}


void JXLImageWriter::save(const char* filename)
{
    JxlEncoderCloseInput(_enc.get());

    std::vector<uint8_t> compressed(64);

    uint8_t* next_out = compressed.data();
    size_t avail_out = compressed.size() - (next_out - compressed.data());
    
    JxlEncoderStatus status = JXL_ENC_NEED_MORE_OUTPUT;

    while (status == JXL_ENC_NEED_MORE_OUTPUT) {
        status = JxlEncoderProcessOutput(_enc.get(), &next_out, &avail_out);

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



// ===========================================================================


JXLImageReader::JXLImageReader(const char* filename)
{
    JxlDecoderStatus status;

    std::vector<uint8_t> jxl_data;

    // Read file
    FILE* file = fopen(filename, "rb");

    fseek(file, 0, SEEK_END);
    const long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    jxl_data.resize(file_size);

    fread(jxl_data.data(), 1, file_size, file);
    fclose(file);

    _dec = JxlDecoderMake(nullptr);
    _runner = JxlThreadParallelRunnerMake(nullptr, JxlThreadParallelRunnerDefaultNumWorkerThreads());

    status = JxlDecoderSetParallelRunner(
        _dec.get(), 
        JxlThreadParallelRunner, 
        _runner.get()
    );

    CHECK_JXL_DEC_STATUS(status);

    status = JxlDecoderSubscribeEvents(
        _dec.get(), 
        JXL_DEC_BASIC_INFO |
        JXL_DEC_BOX |
        // JXL_DEC_COLOR_ENCODING | 
        JXL_DEC_FULL_IMAGE);
    
    CHECK_JXL_DEC_STATUS(status);

    JxlPixelFormat format = {1, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};
    std::vector<uint8_t> box_raw_sgeg;

    JxlDecoderSetInput(_dec.get(), jxl_data.data(), jxl_data.size());

    for (JxlDecoderStatus status_process = JXL_DEC_NEED_MORE_INPUT; 
        status_process != JXL_DEC_FULL_IMAGE; 
        status_process = JxlDecoderProcessInput(_dec.get())) 
    {
        size_t pixel_buffer_sz = 0;
        JxlBoxType box_type;

        switch(status_process) {
            case JXL_DEC_BASIC_INFO:
                status = JxlDecoderGetBasicInfo(_dec.get(), &_basic_info);

                // We need to allocate subframebuffers
                if (_basic_info.num_extra_channels > 0) {
                    _sub_framebuffers.resize(_basic_info.num_extra_channels);
                    _sub_framebuffers_names.resize(_basic_info.num_extra_channels);

                    for (uint32_t i = 0; i < _basic_info.num_extra_channels; i++) {
                        _sub_framebuffers[i].resize(_basic_info.xsize * _basic_info.ysize);

                        JxlExtraChannelInfo extra_info;
                        status = JxlDecoderGetExtraChannelInfo(_dec.get(), i, &extra_info);
                        CHECK_JXL_DEC_STATUS(status);

                        if (extra_info.name_length > 0) {
                            char* channel_name = (char*)malloc(extra_info.name_length + 1);
                            status = JxlDecoderGetExtraChannelName(_dec.get(), i, channel_name, extra_info.name_length + 1);
                            CHECK_JXL_DEC_STATUS(status);

                            _sub_framebuffers_names[i] = channel_name;

                            free(channel_name);
                        }
                    }
                }

                break;
            case JXL_DEC_NEED_IMAGE_OUT_BUFFER:
                std::cout << "JXL_DEC_NEED_IMAGE_OUT_BUFFER" << std::endl;

                // status = JxlDecoderImageOutBufferSize(_dec.get(), &format, &buffer_size);
                // CHECK_JXL_DEC_STATUS(status);

                _main_framebuffer.resize(_basic_info.xsize * _basic_info.ysize);
                pixel_buffer_sz = _basic_info.xsize * _basic_info.ysize * sizeof(float);

                status = JxlDecoderSetImageOutBuffer(
                    _dec.get(), 
                    &format, 
                    _main_framebuffer.data(),
                    pixel_buffer_sz
                );

                CHECK_JXL_DEC_STATUS(status);

                // Do the same for sub framebuffers
                for (uint32_t i = 0; i < _basic_info.num_extra_channels; i++) {
                    status = JxlDecoderSetExtraChannelBuffer(
                        _dec.get(),
                        &format,
                        _sub_framebuffers[i].data(),
                        pixel_buffer_sz,
                        i
                    );

                    CHECK_JXL_DEC_STATUS(status);
                }

                break;
            case JXL_DEC_BOX:
                status = JxlDecoderGetBoxType(_dec.get(), box_type, JXL_TRUE);
                std::cout << "Box detected: " << box_type[0] << box_type[1] << box_type[2] << box_type[3] << std::endl;

                if (box_type[0] == 's' && box_type[1] == 'g' && box_type[2] == 'e' && box_type[3] == 'g') {
                    uint64_t box_size;
                    status = JxlDecoderGetBoxSizeRaw(_dec.get(), &box_size);
                    box_raw_sgeg.resize(box_size);
                    std::cout << "box_size r: " << box_raw_sgeg.size() << std::endl;

                    JxlDecoderSetBoxBuffer(_dec.get(), box_raw_sgeg.data(), box_raw_sgeg.size());
                }
                break;
            case JXL_DEC_NEED_MORE_INPUT:
                std::cerr << "JXL_DEC_NEED_MORE_INPUT" << std::endl;
                break;
            case JXL_DEC_ERROR:
                std::cerr << "JXL_DEC_ERROR!" << std::endl;
                exit(1);
                break;
            default:
                std::cout << "Unkown decoder status: " << status_process << std::endl;
                break;
        }
    }

    JxlDecoderReleaseInput(_dec.get());
    
    _sgeg_box = SGEG_box(box_raw_sgeg);

    CHECK_JXL_DEC_STATUS(status);
}


JXLImageReader::~JXLImageReader() {
    // JxlDecoderDestroy(_dec.get());
}


uint32_t JXLImageReader::width() const
{
    return _basic_info.xsize;
}


uint32_t JXLImageReader::height() const
{
    return _basic_info.ysize;
}


uint32_t JXLImageReader::n_subframebuffers() const
{
    return _sub_framebuffers.size();
}


SGEG_box JXLImageReader::get_sgeg() const
{
    return _sgeg_box;
}


void JXLImageReader::getMainFramebuffer(std::vector<float>& framebuffer) const
{
    framebuffer.resize(_main_framebuffer.size());

    memcpy(framebuffer.data(), _main_framebuffer.data(), sizeof(float) * width() * height());
}


void JXLImageReader::getSubFramebuffer(std::vector<float>& framebuffer, size_t index) const
{
    if (index >= _sub_framebuffers.size()) {
        std::cerr << "Requested index out of range" << std::endl;
        return;
    }
    
    framebuffer.resize(_sub_framebuffers[index].size());

    memcpy(framebuffer.data(), _sub_framebuffers[index].data(), sizeof(float) * width() * height());
}


void JXLImageReader::print_basic_info() const
{
    std::cout << "┌                   width: " << _basic_info.xsize << std::endl;
    std::cout << "├                  height: " << _basic_info.ysize << std::endl;
    std::cout << "├         bits per sample: " << _basic_info.bits_per_sample << std::endl;
    std::cout << "├exponent bits per sample: " << _basic_info.exponent_bits_per_sample << std::endl;
    std::cout << "├        intensity target: " << _basic_info.intensity_target << std::endl;
    std::cout << "├ relative to max display: " << _basic_info.relative_to_max_display << std::endl;
    std::cout << "├            linear below: " << _basic_info.linear_below << std::endl;
    std::cout << "├   uses original profile: " << _basic_info.uses_original_profile << std::endl;
    std::cout << "├            have preview: " << _basic_info.have_preview << std::endl;
    std::cout << "├          have animation: " << _basic_info.have_animation << std::endl;
    std::cout << "├             orientation: " << _basic_info.orientation << std::endl;
    std::cout << "├      num color channels: " << _basic_info.num_color_channels << std::endl;
    std::cout << "├      num extra channels: " << _basic_info.num_extra_channels << std::endl;
    std::cout << "├              alpha bits: " << _basic_info.alpha_bits << std::endl;
    std::cout << "├     alpha exponent bits: " << _basic_info.alpha_exponent_bits << std::endl;
    std::cout << "└     alpha premultiplied: " << _basic_info.alpha_premultiplied << std::endl;


    for (uint32_t i = 0; i < _basic_info.num_extra_channels; i++) {
        std::cout << "\t┌ Channel #" << i << std::endl;
        std::cout << "\t└  name: " << _sub_framebuffers_names[i] << std::endl; 
    }

    for (uint32_t i = 0; i < _sgeg_box.n_moments; i++) {
        std::cout << "xmin: " << _sgeg_box.moment_min[i] << std::endl;
        std::cout << "xmax: " << _sgeg_box.moment_max[i] << std::endl;
    }
}
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

#include "JXLImage.h"
#include "Util.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <sstream>

#include <half.h>

// ----------------------------------------------------------------------------

#define JXL_MAX_FRAMEBUFFERS 1
//256

#define CHECK_JXL_ENC_STATUS(status)                                          \
    if (JXL_ENC_SUCCESS != (status)) {                                        \
        std::stringstream err_msg;                                            \
        err_msg << "[ERR] JXL_ENCODER: "                                      \
                << __FILE__ << ":" << __LINE__                                \
                << "> " << status;                                            \
        throw std::runtime_error(err_msg.str());                              \
    }                                                                         \


#define CHECK_JXL_DEC_STATUS(status)                                          \
    if (JXL_DEC_SUCCESS != (status)) {                                        \
        std::stringstream err_msg;                                            \
        err_msg << "[ERR] JXL_DECODER: "                                      \
                << __FILE__ << ":" << __LINE__                                \
                << "> " << status;                                            \
        throw std::runtime_error(err_msg.str());                              \
    }                                                                         \

// ----------------------------------------------------------------------------

JXLFramebuffer::JXLFramebuffer(
    uint32_t width, uint32_t height,
    uint32_t n_color_channels,
    std::pair<int, int> n_bits_per_sample,
    uint32_t subsampling_factor,
    float    framedistance,
    const char* name)
    : _name(nullptr)
    , _n_bits_per_sample(n_bits_per_sample)
    , _subsampling_factor(subsampling_factor)
    , _framedistance(framedistance)
    , _pixel_format(
        {
            /* .num_channels = */ n_color_channels,
            /* .data_type    = */ JXL_TYPE_FLOAT,
            /* .endianness   = */ JXL_NATIVE_ENDIAN,
            /* .align        = */ 0
        })
    , _pixel_data(width * height * n_color_channels)
{
    assert(_pixel_data.size() > 0);

    if (name != nullptr) { setName(name); }
}



JXLFramebuffer::JXLFramebuffer(
    const std::vector<float>& framebuffer,
    uint32_t n_color_channels,
    std::pair<int, int> n_bits_per_sample,
    uint32_t subsampling_factor,
    float    framedistance,
    const char* name)
    : _name(nullptr)
    , _n_bits_per_sample(n_bits_per_sample)
    , _subsampling_factor(subsampling_factor)
    , _framedistance(framedistance)
    , _pixel_format(
        {
            /* .num_channels = */ n_color_channels,
            /* .data_type    = */ JXL_TYPE_FLOAT,
            /* .endianness   = */ JXL_NATIVE_ENDIAN,
            /* .align        = */ 0
        })
    , _pixel_data(framebuffer)
{
    assert(_pixel_data.size() > 0);

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


void JXLFramebuffer::dump(FILE* stream) const
{
    const uint32_t len_name = (_name) ? std::strlen(_name) : 0;
    std::fwrite(&len_name, sizeof(uint32_t), 1, stream);

    if (len_name > 0) {
        std::fwrite(_name, sizeof(char), len_name, stream);
    }

    std::fwrite(&_n_bits_per_sample.first, sizeof(uint32_t), 1, stream);
    std::fwrite(&_n_bits_per_sample.second, sizeof(uint32_t), 1, stream);
    std::fwrite(&_subsampling_factor, sizeof(uint32_t), 1, stream);
    std::fwrite(&_pixel_format, sizeof(JxlPixelFormat), 1, stream);

    const uint32_t n_pixels = _pixel_data.size();
    assert(n_pixels > 0);

    std::fwrite(&n_pixels, sizeof(uint32_t), 1, stream);
    std::fwrite(_pixel_data.data(), sizeof(float), n_pixels, stream);
}


JXLFramebuffer* JXLFramebuffer::read_dump(FILE* stream)
{
    JXLFramebuffer* fb = new JXLFramebuffer(1, 1);

    uint32_t len_name;
    std::fread(&len_name, sizeof(uint32_t), 1, stream);

    if (len_name > 0) {
        fb->_name = new char[len_name + 1];
        std::fread(fb->_name, sizeof(char), len_name, stream);
    } else {
        fb->_name = nullptr;
    }

    std::fread(&(fb->_n_bits_per_sample.first), sizeof(uint32_t), 1, stream);
    std::fread(&(fb->_n_bits_per_sample.second), sizeof(uint32_t), 1, stream);
    std::fread(&(fb->_subsampling_factor), sizeof(uint32_t), 1, stream);
    std::fread(&(fb->_pixel_format), sizeof(JxlPixelFormat), 1, stream);

    uint32_t n_pixels;

    std::fread(&n_pixels, sizeof(uint32_t), 1, stream);
    fb->_pixel_data.resize(n_pixels);
    std::fread(fb->_pixel_data.data(), sizeof(float), n_pixels, stream);

    return fb;
}

// ===========================================================================



JXLImage::JXLImage(const char* filename)
    : _n_parts(1) // At start, we assume there is only a single file
{
    load(filename);
}


JXLImage::JXLImage(const std::string& filename)
    : _n_parts(1) // At start, we assume there is only a single file
{
    load(filename.c_str());
}


JXLImage::JXLImage(uint32_t width, uint32_t height)
    : _width(width)
    , _height(height)
    , _n_parts(1)
{}


JXLImage::~JXLImage()
{
    for (size_t i = 0; i < _framebuffers.size(); i++) {
        delete _framebuffers[i];
    }
}


void JXLImage::setBox(const SGEGBox& box)
{
    _sgeg_box = box;
    _sgeg_box.n_parts = _n_parts;
}


size_t JXLImage::appendFramebuffer(
    const std::vector<float>& framebuffer,
    uint32_t n_channels,
    std::pair<int, int> enc_bits_per_sample,
    uint32_t enc_subsampling_factor,
    float    enc_framedistance,
    const char* name)
{
    assert(framebuffer.size() == n_channels * _width * _height);

    JXLFramebuffer *fb = new JXLFramebuffer(
        framebuffer,
        n_channels,
        enc_bits_per_sample,
        enc_subsampling_factor,
        enc_framedistance,
        name
    );

    _framebuffers.push_back(fb);

    _n_parts = std::ceil((float)_framebuffers.size() / (float)JXL_MAX_FRAMEBUFFERS);

    return _framebuffers.size() - 1;
}


void JXLImage::write(const char* filename, int effort) const {
    // Currently, JXL supports up to 256 framebuffers per image
    // we may need more than that so, in such scenario, we are
    // going to write multiple images
    std::string base, ext;
    Util::split_extension(filename, base, ext);

    for (uint32_t part = 0; part < _n_parts; part++) {
        std::string curr_filename;

        if (part == 0) {
            curr_filename = base + ext;
        } else {
            std::stringstream ss;
            ss << base << "_" << part << ext;
            curr_filename = ss.str();
        }

        std::vector<float> quantized_framebuffer;

        const int start_framebuffer_idx = part * JXL_MAX_FRAMEBUFFERS;
        const int end_framebuffer_idx = std::min(start_framebuffer_idx + JXL_MAX_FRAMEBUFFERS, (int)_framebuffers.size()) - 1;
        const int n_framebuffers = end_framebuffer_idx - start_framebuffer_idx + 1;

        assert(start_framebuffer_idx >= 0);
        assert(end_framebuffer_idx >= start_framebuffer_idx);
        assert(end_framebuffer_idx < (int)_framebuffers.size());
        assert(n_framebuffers <= JXL_MAX_FRAMEBUFFERS);

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

        // --------------------------------------------------------------------

        // We save the box just for the 1st image
        if (part == 0) {
            status = JxlEncoderUseBoxes(enc.get());
            CHECK_JXL_ENC_STATUS(status);

            char tp[4] = {'s','g','e','g'};

            std::vector<uint8_t> raw_box;
            _sgeg_box.getRaw(raw_box);

            status = JxlEncoderAddBox(enc.get(), tp, raw_box.data(), raw_box.size(), JXL_FALSE);
            CHECK_JXL_ENC_STATUS(status);
        }

        // ====================================================================
        // Main framebuffer
        // ====================================================================

        JxlEncoderInitBasicInfo(&basic_info);

        basic_info.xsize                    = _width;
        basic_info.ysize                    = _height;
        basic_info.num_extra_channels       = n_framebuffers - 1;
        basic_info.num_color_channels       = _framebuffers[start_framebuffer_idx]->getNColorChannels();
        basic_info.bits_per_sample          = _framebuffers[start_framebuffer_idx]->getBitsPerSample();
        basic_info.exponent_bits_per_sample = _framebuffers[start_framebuffer_idx]->getExponentBitsPerSample();
        basic_info.uses_original_profile    = JXL_TRUE;

        status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);
        CHECK_JXL_ENC_STATUS(status);

        // --------------------------------------------------------------------

        JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

        // Set compression quality
        JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, effort);
        CHECK_JXL_ENC_STATUS(status);

        const bool encodes_lossless = _framebuffers[start_framebuffer_idx]->getFramedistance() == 0;

        if (encodes_lossless) {
            status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
        } else {
            status = JxlEncoderSetFrameLossless(frame_settings, JXL_FALSE);
            status = JxlEncoderSetFrameDistance(frame_settings, _framebuffers[start_framebuffer_idx]->getFramedistance());
        }
        CHECK_JXL_ENC_STATUS(status);

        status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_RESAMPLING, _framebuffers[start_framebuffer_idx]->getDownsamplingFactor());
        CHECK_JXL_ENC_STATUS(status);

        // --------------------------------------------------------------------

        switch (_framebuffers[start_framebuffer_idx]->getNColorChannels()) {
            case 1:
                color_encoding.color_space = JXL_COLOR_SPACE_GRAY;
                break;
            case 3:
            case 4:
                color_encoding.color_space = JXL_COLOR_SPACE_RGB;
                break;
            default:
                throw std::runtime_error("Unknown color space matching number of channels");
                break;
        }

        color_encoding.white_point       = JXL_WHITE_POINT_D65;
        color_encoding.primaries         = JXL_PRIMARIES_SRGB;
        color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
        color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;

        // TODO set name

        status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);
        CHECK_JXL_ENC_STATUS(status);

        const JxlPixelFormat format = _framebuffers[start_framebuffer_idx]->getPixelFormat();
        const size_t data_size      = _framebuffers[start_framebuffer_idx]->getPixelData().size() * sizeof(float);

        // FIXME: There is a bug in the current JXL implementation: it ignores
        //        the quantization when the compression is not lossless.
        //        This issue is addressed with a hack.
        if (encodes_lossless || basic_info.exponent_bits_per_sample != 0) {
            // When we are dealing with lossless file and halfs, libjxl
            // triggers an error if the half rounding cause imprecision.
            // So we explicitely cast to half and back to float before
            // providing the buffer.
            if (encodes_lossless && basic_info.exponent_bits_per_sample == 5
             && basic_info.bits_per_sample == 16) {
                quantized_framebuffer.resize(_width * _height);

                #pragma omp parallel for
                for (size_t i = 0; i < _width * _height; i++) {
                    quantized_framebuffer[i] =
                        imath_half_to_float(
                            imath_float_to_half(_framebuffers[start_framebuffer_idx]->getPixelData()[i]
                        )
                    );
                }

                status = JxlEncoderAddImageFrame(
                  frame_settings,
                  &format,
                  quantized_framebuffer.data(),
                  data_size
                );
            } else {
                status = JxlEncoderAddImageFrame(
                    frame_settings,
                    &format,
                    _framebuffers[start_framebuffer_idx]->getPixelData().data(),
                    data_size
                );
            }
        } else {
            quantized_framebuffer.resize(_width * _height);

            #pragma omp parallel for
            for (size_t i = 0; i < _width * _height; i++) {
                quantized_framebuffer[i] = Util::quantize_dequantize(
                    _framebuffers[start_framebuffer_idx]->getPixelData()[i],
                    basic_info.bits_per_sample
                );
            }

            status = JxlEncoderAddImageFrame(
                frame_settings,
                &format,
                quantized_framebuffer.data(),
                data_size
            );
        }

        CHECK_JXL_ENC_STATUS(status);

        // ====================================================================
        // Subframebuffers
        // ====================================================================

        for (int i = 1; i < n_framebuffers; i++) {
            const JXLFramebuffer* fb    = _framebuffers[start_framebuffer_idx + i];
            const JxlPixelFormat format = fb->getPixelFormat();
            const size_t data_size      = fb->getPixelDataConst().size() * sizeof(float);

            JxlExtraChannelInfo extra_info;
            JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);

            extra_info.bits_per_sample          = fb->getBitsPerSample();
            extra_info.exponent_bits_per_sample = fb->getExponentBitsPerSample();

            JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

            status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, effort);
            CHECK_JXL_ENC_STATUS(status);

            // Set compression quality
            // FIXME: Not supported yet in the version of libjxl used (0.7)
            //        at the time of writing this code
            // JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EXTRA_CHANNEL_RESAMPLING, subsampling);
            // status = JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EXTRA_CHANNEL_RESAMPLING, 4);
            // CHECK_JXL_ENC_STATUS(status);

            // FIXME: This does absolutely nothing at the time this code was
            //        written. In a near future, we expect a support for custom
            //        compression ratio per subimage in libjxl.
            //        Also, when implementing, worth checking if the
            //        quantization works as expected.
            if (fb->getFramedistance() > 0) {
                status = JxlEncoderSetFrameLossless(frame_settings, JXL_FALSE);
                status = JxlEncoderSetFrameDistance(frame_settings, fb->getFramedistance());
            } else {
                status = JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE);
            }
            CHECK_JXL_ENC_STATUS(status);

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
        size_t avail_out  = compressed.size();

        status = JXL_ENC_NEED_MORE_OUTPUT;

        while (status == JXL_ENC_NEED_MORE_OUTPUT) {
            status = JxlEncoderProcessOutput(enc.get(), &next_out, &avail_out);

            if (status == JXL_ENC_NEED_MORE_OUTPUT) {
                size_t offset = next_out - compressed.data();
                compressed.resize(compressed.size() * 2);

                next_out  = compressed.data() + offset;
                avail_out = compressed.size() - offset;
            }
        }

        CHECK_JXL_ENC_STATUS(status);

        compressed.resize(next_out - compressed.data());

        // Write file
        std::FILE* file = std::fopen(curr_filename.data(), "wb");

        if (!file) {
            throw std::runtime_error("Could not open file " + curr_filename + " for writing.");
        }

        std::fwrite(compressed.data(), sizeof(uint8_t), compressed.size(), file);
        std::fclose(file);
    }
}


void JXLImage::write(const std::string& filename, int effort) const
{
    write(filename.c_str(), effort);
}


void JXLImage::dump(const char* filename) const
{
    std::FILE* f = std::fopen(filename, "wb");

    if (!f) {
        throw std::runtime_error("Could not create an JXLImage dump file");
        return;
    }

    std::fwrite(&_width, sizeof(uint32_t), 1, f);
    std::fwrite(&_height, sizeof(uint32_t), 1, f);

    std::vector<uint8_t> sgeg_data;
    _sgeg_box.getRaw(sgeg_data);
    const uint32_t sgeg_data_size = sgeg_data.size();
    std::fwrite(&sgeg_data_size, sizeof(uint32_t), 1, f);
    std::fwrite(sgeg_data.data(), sizeof(uint8_t), sgeg_data_size, f);

    const uint32_t n_framebuffers = _framebuffers.size();
    std::fwrite(&n_framebuffers, sizeof(uint32_t), 1, f);

    for (const JXLFramebuffer* fb: _framebuffers) {
        fb->dump(f);
    }

    std::fclose(f);
}


void JXLImage::dump(const std::string& filename) const
{
    dump(filename.c_str());
}


JXLImage* JXLImage::read_dump(const char* filename)
{
    std::FILE *f = std::fopen(filename, "rb");

    if (!f) {
        throw std::runtime_error("Could not read JXLImage dump file");
        return nullptr;
    }

    JXLImage* image = new JXLImage(1, 1);

    std::fread(&(image->_width), sizeof(uint32_t), 1, f);
    std::fread(&(image->_height), sizeof(uint32_t), 1, f);

    uint32_t sgeg_data_size;
    std::fread(&sgeg_data_size, sizeof(uint32_t), 1, f);
    std::vector<uint8_t> sgeg_data(sgeg_data_size);
    std::fread(sgeg_data.data(), sizeof(uint8_t), sgeg_data_size, f);

    image->_sgeg_box = SGEGBox(sgeg_data);

    uint32_t n_framebuffers;

    std::fread(&n_framebuffers, sizeof(uint32_t), 1, f);

    image->_framebuffers.reserve(n_framebuffers);

    for (size_t i = 0; i < n_framebuffers; i++) {
        image->_framebuffers.push_back(JXLFramebuffer::read_dump(f));
    }

    std::fclose(f);

    return image;
}


void JXLImage::load(const char* filename)
{
    // Currently, JXL supports up to 256 framebuffers per image
    // we may need have need more when writing the image...
    // So extra images must be taken care of.

    std::string base, ext;
    Util::split_extension(filename, base, ext);

    std::vector<uint8_t> box_raw_sgeg;

    // At start, we assume there is only a single file, the SGEG
    // box will tell us uppon decompression of the main file if
    // additional files exist

    for (uint32_t part = 0; part < _n_parts; part++) {
        std::string curr_filename;

        if (part == 0) {
            curr_filename = base + ext;
        } else {
            std::stringstream ss;
            ss << base << "_" << part << ext;
            curr_filename = ss.str();
        }

        const uint32_t start_framebuffer_idx = part * JXL_MAX_FRAMEBUFFERS;
        uint32_t n_framebuffer_part = 0;

        JxlDecoderStatus status;

        std::vector<uint8_t> jxl_data;

        JxlDecoderPtr              dec;
        JxlThreadParallelRunnerPtr runner;

        // Read file
        std::FILE* file = fopen(curr_filename.data(), "rb");

        if (file == NULL) {
            throw std::runtime_error("Could not open " + curr_filename);
        }

        std::fseek(file, 0, SEEK_END);
        const long file_size = ftell(file);
        std::fseek(file, 0, SEEK_SET);

        jxl_data.resize(file_size);

        std::fread(jxl_data.data(), 1, file_size, file);
        std::fclose(file);

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
                            std::make_pair(basic_info.bits_per_sample, basic_info.exponent_bits_per_sample)
                        ));

                        // Extra layer metadata
                        n_framebuffer_part = basic_info.num_extra_channels + 1;

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
                                std::make_pair(extra_info.bits_per_sample, extra_info.exponent_bits_per_sample),
                                1,
                                0, // TODO: The current API does not gives the used framedistance
                                (extra_info.name_length > 0) ? layer_name.data() : nullptr
                            ));
                        }
                    }
                    break;

                case JXL_DEC_NEED_IMAGE_OUT_BUFFER:
                    {
                        JxlPixelFormat format = _framebuffers[start_framebuffer_idx]->getPixelFormat();

                        status = JxlDecoderSetImageOutBuffer(
                            dec.get(),
                            &format,
                            _framebuffers[start_framebuffer_idx]->getPixelData().data(),
                            _framebuffers[start_framebuffer_idx]->getSizeBytes()
                        );
                        CHECK_JXL_DEC_STATUS(status);

                        for (size_t i = 1; i < n_framebuffer_part; i++) {
                            format = _framebuffers[start_framebuffer_idx + i]->getPixelFormat();

                            status = JxlDecoderSetExtraChannelBuffer(
                                dec.get(),
                                &format,
                                _framebuffers[start_framebuffer_idx + i]->getPixelData().data(),
                                _framebuffers[start_framebuffer_idx + i]->getSizeBytes(),
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
                            // The SGEG box must be present only on part 0
                            assert(part == 0);

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
                    std::stringstream err_msg;
                    err_msg << "Unknown decoder status: " << status_process;
                    throw std::runtime_error(err_msg.str());
                    break;
            }
        }

        JxlDecoderReleaseInput(dec.get());
        CHECK_JXL_DEC_STATUS(status);

        // For part 0, ensure we have a SGEG box
        if (part == 0) {
            assert(box_raw_sgeg.size() > 0);
            _sgeg_box = SGEGBox(box_raw_sgeg);
            _n_parts = _sgeg_box.n_parts;
        }
    }
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

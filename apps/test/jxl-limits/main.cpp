#include <iostream>
#include <vector>

#include <jxl/encode_cxx.h>
#include <jxl/thread_parallel_runner_cxx.h>

// ----------------------------------------------------------------------------

#define CHECK_JXL_ENC_STATUS(status)                                          \
    if (JXL_ENC_SUCCESS != status) {                                          \
        std::cerr << "[ERR] At: "                                             \
                  << __FILE__  << ":" << __LINE__                             \
                  << " " << status << std::endl;                              \
    }                                                                         \

// ----------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    size_t width = 128;
    size_t height = 64;
    size_t n_extra_channels = 256;
    uint32_t n_color_channels = 3;

    std::vector<uint8_t> main_image_data(width * height * n_color_channels);
    std::vector<uint8_t> sub_image_data(width * height);

    const size_t main_data_size = main_image_data.size() * sizeof(uint8_t);
    const size_t sub_data_size  = sub_image_data.size() * sizeof(uint8_t);

    const JxlPixelFormat main_format = {n_color_channels, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};
    const JxlPixelFormat sub_format  = {1, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};

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

    JxlEncoderInitBasicInfo(&basic_info);

    basic_info.xsize                    = width;
    basic_info.ysize                    = height;
    basic_info.num_extra_channels       = n_extra_channels;
    basic_info.num_color_channels       = n_color_channels;
    basic_info.bits_per_sample          = 8;
    basic_info.exponent_bits_per_sample = 0;
    basic_info.uses_original_profile    = JXL_TRUE;

    status = JxlEncoderSetBasicInfo(enc.get(), &basic_info);

    CHECK_JXL_ENC_STATUS(status);

    color_encoding.color_space       = (n_color_channels == 1) ? JXL_COLOR_SPACE_GRAY: JXL_COLOR_SPACE_RGB;
    color_encoding.white_point       = JXL_WHITE_POINT_D65;
    color_encoding.primaries         = JXL_PRIMARIES_SRGB;
    color_encoding.transfer_function = JXL_TRANSFER_FUNCTION_LINEAR;
    color_encoding.rendering_intent  = JXL_RENDERING_INTENT_PERCEPTUAL;
    
    status = JxlEncoderSetColorEncoding(enc.get(), &color_encoding);

    CHECK_JXL_ENC_STATUS(status);

    JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

    status = JxlEncoderAddImageFrame(
        frame_settings, 
        &main_format,
        main_image_data.data(),
        main_data_size
    );

    CHECK_JXL_ENC_STATUS(status);

    for (size_t i = 0; i < n_extra_channels; i++) {
        JxlExtraChannelInfo extra_info;
        JxlEncoderInitExtraChannelInfo(JXL_CHANNEL_OPTIONAL, &extra_info);

        extra_info.bits_per_sample          = 8;
        extra_info.exponent_bits_per_sample = 0;

        JxlEncoderFrameSettings* frame_settings = JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

        status = JxlEncoderSetExtraChannelInfo(enc.get(), i, &extra_info);
        CHECK_JXL_ENC_STATUS(status);

        status = JxlEncoderSetExtraChannelBuffer(
            frame_settings, 
            &sub_format, 
            sub_image_data.data(),
            sub_data_size,
            i);

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
    FILE* file = fopen("test.jxl", "wb");
    fwrite(compressed.data(), sizeof(uint8_t), compressed.size(), file);
    fclose(file);


    return 0;
}
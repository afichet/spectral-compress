set(PUBLIC_HEADERS
    moments.h
    moments_image.h
    moments_error.h
    curve_quantization.h
    curve_compression.h
    JXLImage.h
    EXRImage.h
    EXRSpectralImage.h
    EXRArrayStream.h
    SGEGBox.h
    SpectrumConverter.h
    Util.h
    spectral_compression_type.h
)

add_library(common STATIC
    moments.cpp
    moments_image.cpp
    moments_error.cpp
    curve_quantization.cpp
    curve_compression.cpp
    JXLImage.cpp
    EXRImage.cpp
    EXRSpectralImage.cpp
    EXRArrayStream.cpp
    SGEGBox.cpp
    SpectrumConverter.cpp
    Util.cpp
    EXRUtils.cpp
)

set_target_properties(common PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_include_directories(common PUBLIC .)

find_package(OpenMP)

if (OpenMP_FOUND OR OpenMP_CXX_FOUND)
   target_link_libraries(common PUBLIC OpenMP::OpenMP_CXX)
endif()

if (MSVC)
    target_compile_options(common PUBLIC /W3)
else()
    target_link_libraries(common PUBLIC m)
    target_compile_options(common PUBLIC -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(common PUBLIC jxl jxl_threads OpenEXR::OpenEXR Imath::Imath Eigen3::Eigen)

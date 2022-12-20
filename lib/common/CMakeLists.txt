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
)

set_target_properties(common PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_include_directories(common PUBLIC .)

find_package(OpenMP)

if (OpenMP_FOUND OR OpenMP_CXX_FOUND)
   target_link_libraries(common PUBLIC OpenMP::OpenMP_CXX)
endif()

# OpenEXR 3.0 provides its own find_package
find_package(OpenEXR CONFIG)
find_package(Imath CONFIG)

if (OpenEXR_FOUND)
    set(OpenEXR_LIBRARIES OpenEXR::OpenEXR Imath::Imath)
else()
    # For OpenEXR 2.0 we rely on our own find_package
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(OpenEXR REQUIRED)
endif()

if (MSVC)
    target_compile_options(common PUBLIC /W3)
else()
    target_link_libraries(common PUBLIC m)
    target_compile_options(common PUBLIC -Wall -Wextra -Wpedantic)
endif()

target_include_directories(common PUBLIC ${OpenEXR_INCLUDE_DIR})
target_link_libraries(common PUBLIC ${OpenEXR_LIBRARIES})
target_link_libraries(common PUBLIC jxl jxl_threads Eigen3::Eigen)

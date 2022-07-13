add_executable(decompress main.cpp)

include_directories("../../3rdparty/heif/libheif")

target_link_libraries(decompress PRIVATE heif EXRSpectralImage moments)
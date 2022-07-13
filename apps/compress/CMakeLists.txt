add_executable(compress main.cpp)

include_directories("../../3rdparty/heif/libheif")


target_link_libraries(compress PRIVATE heif EXRSpectralImage moments)
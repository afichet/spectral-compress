add_executable(exr-peek-spectrum main.cpp)

target_link_libraries(exr-peek-spectrum PRIVATE common lodepng)
target_include_directories(exr-peek-spectrum PRIVATE ../../../3rdparty/tclap/include)

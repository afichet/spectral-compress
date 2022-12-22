add_executable(exr-diff main.cpp)

target_link_libraries(exr-diff PRIVATE common lodepng)
target_include_directories(exr-diff PRIVATE ../../../3rdparty/tclap/include)

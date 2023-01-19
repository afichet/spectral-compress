add_executable(exr-strip-rgb main.cpp)

target_link_libraries(exr-strip-rgb PRIVATE common lodepng)
target_include_directories(exr-strip-rgb PRIVATE ../../../3rdparty/tclap/include)

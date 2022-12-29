add_executable(exr-png main.cpp)

target_link_libraries(exr-png PRIVATE common lodepng)
target_include_directories(exr-png PRIVATE ../../../3rdparty/tclap/include)

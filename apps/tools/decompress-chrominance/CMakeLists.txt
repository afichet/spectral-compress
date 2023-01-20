add_executable(decompress-chrominance main.cpp)

target_link_libraries(decompress-chrominance PRIVATE common)
target_include_directories(decompress-chrominance PRIVATE ../../../3rdparty/tclap/include)

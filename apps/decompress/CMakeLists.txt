add_executable(decompress main.cpp)

target_link_libraries(decompress PRIVATE common)
target_include_directories(decompress PRIVATE ../../3rdparty/tclap/include)

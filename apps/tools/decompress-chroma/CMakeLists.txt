add_executable(decompress-chroma main.cpp)

target_link_libraries(decompress-chroma PRIVATE common)
target_include_directories(decompress-chroma PRIVATE ../../../3rdparty/tclap/include)

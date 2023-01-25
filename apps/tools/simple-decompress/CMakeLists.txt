add_executable(simple-decompress main.cpp)

target_link_libraries(simple-decompress PRIVATE common)
target_include_directories(simple-decompress PRIVATE ../../../3rdparty/tclap/include)

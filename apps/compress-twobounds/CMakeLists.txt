add_executable(compress-twobounds main.cpp)

target_link_libraries(compress-twobounds PRIVATE common)
target_include_directories(compress-twobounds PRIVATE ../../3rdparty/tclap/include)

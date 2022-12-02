add_executable(compress-upperbound main.cpp)

target_link_libraries(compress-upperbound PRIVATE common)
target_include_directories(compress-upperbound PRIVATE ../../3rdparty/tclap/include)

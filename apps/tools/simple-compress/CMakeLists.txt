add_executable(simple-compress main.cpp)
target_link_libraries(simple-compress PRIVATE common)
target_include_directories(simple-compress PRIVATE ../../../3rdparty/tclap/include)

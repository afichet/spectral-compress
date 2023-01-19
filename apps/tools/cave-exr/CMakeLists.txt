add_executable(cave-exr main.cpp)

target_link_libraries(cave-exr PRIVATE common png)
target_include_directories(cave-exr PRIVATE ../../../3rdparty/tclap/include)

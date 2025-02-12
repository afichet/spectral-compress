add_executable(cave-exr main.cpp)

find_package(PNG)

target_link_libraries(cave-exr PRIVATE common PNG::PNG)
target_include_directories(cave-exr PRIVATE ../../../3rdparty/tclap/include)

add_executable(exr-change-compression main.cpp)

target_link_libraries(exr-change-compression PRIVATE common)
target_include_directories(exr-change-compression PRIVATE ../../../3rdparty/tclap/include)

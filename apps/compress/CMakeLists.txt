add_executable(compress main.cpp stats.cpp cmd_constraints.cpp)

target_link_libraries(compress PRIVATE common)
target_include_directories(compress PRIVATE ../../3rdparty/tclap/include)

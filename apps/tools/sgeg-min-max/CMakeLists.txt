add_executable(sgeg-min-max main.cpp)
target_link_libraries(sgeg-min-max PRIVATE common)
target_include_directories(sgeg-min-max PRIVATE ../../../3rdparty/tclap/include)

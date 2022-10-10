add_library(lodepng lodepng.cpp)
set_target_properties(lodepng PROPERTIES PUBLIC_HEADER lodepng.h)
target_include_directories(lodepng PUBLIC .)


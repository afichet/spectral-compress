set(PUBLIC_HEADERS
    moments.h
    moments_image.h
    JXLImage.h
    SGEG_box.h
)

add_library(moments STATIC 
    moments.cpp
    moments_image.cpp
    JXLImage.cpp
    SGEG_box.cpp
)

set_target_properties(moments PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_include_directories(moments PUBLIC .)

find_package(OpenMP)

if (OpenMP_FOUND OR OpenMP_CXX_FOUND)
   target_link_libraries(moments PUBLIC OpenMP::OpenMP_CXX)
endif()

if (MSVC)
    target_compile_options(moments PUBLIC /W3)
else()
    target_link_libraries(moments PUBLIC m)
    target_compile_options(moments PUBLIC -Wall -Wextra -Wpedantic)
endif()


target_link_libraries(moments PUBLIC jxl jxl_threads)
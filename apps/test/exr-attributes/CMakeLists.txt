# OpenEXR 3.0 provides its own find_package
find_package(OpenEXR CONFIG)
find_package(Imath CONFIG)

if (OpenEXR_FOUND)
    set(OpenEXR_LIBRARIES OpenEXR::OpenEXR Imath::Imath)
else()
    # For OpenEXR 2.0 we rely on our own find_package
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(OpenEXR REQUIRED)
endif()

add_executable(exr-attributes main.cpp)

target_include_directories(exr-attributes PRIVATE ${OpenEXR_INCLUDE_DIR})
target_link_libraries(exr-attributes PRIVATE ${OpenEXR_LIBRARIES})
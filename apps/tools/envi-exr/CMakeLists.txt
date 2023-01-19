add_executable(envi-exr
    main.cpp
    external/ENVISpectralImage.cpp
)

target_include_directories(envi-exr PRIVATE external)

target_link_libraries(envi-exr PRIVATE common)

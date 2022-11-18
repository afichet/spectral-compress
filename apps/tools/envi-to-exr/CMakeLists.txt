add_executable(envi-to-exr 
    main.cpp
    external/ENVISpectralImage.cpp
)

target_include_directories(envi-to-exr PRIVATE external)

target_link_libraries(envi-to-exr PRIVATE common)
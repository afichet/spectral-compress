add_executable(export-layers
    main.cpp
    lodepng.cpp
    )

target_link_libraries(export-layers PRIVATE EXRSpectralImage common)
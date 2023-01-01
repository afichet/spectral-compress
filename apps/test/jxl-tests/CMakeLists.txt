add_executable(jxl-tests main.cpp)

# target_link_libraries(jxl-foo PRIVATE jxl jxl_threads)
target_link_libraries(jxl-tests PRIVATE common)

# find_package(PkgConfig)
# pkg_check_modules(Jxl REQUIRED IMPORTED_TARGET libjxl)
# pkg_check_modules(JxlThreads REQUIRED IMPORTED_TARGET libjxl_threads)

# target_link_libraries(jxl-foo PkgConfig::Jxl PkgConfig::JxlThreads)

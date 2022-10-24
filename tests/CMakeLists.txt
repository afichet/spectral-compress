add_executable(moment_test moment_test.cpp)
target_link_libraries(moment_test GTest::gtest_main common)

add_executable(util_test util_test.cpp)
target_link_libraries(util_test GTest::gtest_main GTest::gmock common)

add_executable(spectrumconverter_test spectrumconverter_test.cpp)
target_link_libraries(spectrumconverter_test GTest::gtest_main common)

include(GoogleTest)

gtest_discover_tests(moment_test)
gtest_discover_tests(util_test)
gtest_discover_tests(spectrumconverter_test)
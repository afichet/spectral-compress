cmake_minimum_required(VERSION 3.1.1)

project(spectral-compress)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(3rdparty)

set(ENABLE_TEST YES)

# Add flags to the compiler in case we want coverage report
if (ENABLE_TEST)
    if(UNIX)
        set(CMAKE_CXX_FLAGS           ${CMAKE_CXX_FLAGS}           "-fprofile-arcs -ftest-coverage")
        set(CMAKE_C_FLAGS             ${CMAKE_C_FLAGS}             "-fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS    ${CMAKE_EXE_LINKER_FLAGS}    "-fprofile-arcs -ftest-coverage")
        set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "-fprofile-arcs -ftest-coverage")
    endif()  
endif()

add_subdirectory(lib)
add_subdirectory(apps)

# Add Google Test package and test directory
if (ENABLE_TEST)
    message("-- Including Tests")
    include(FetchContent)
        FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
        )

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)

    enable_testing()

    add_subdirectory(tests)
endif()
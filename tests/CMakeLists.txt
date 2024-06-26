cmake_minimum_required(VERSION 3.14)

project(test_magneto LANGUAGES CXX)

# A hack to disable warnings for any targets except magneto
if(DEFINED CMAKE_CXX_CLANG_TIDY)
    set(CMAKE_CXX_CLANG_TIDY_save "${CMAKE_CXX_CLANG_TIDY}")
    set(CMAKE_CXX_CLANG_TIDY "")
endif()

include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY "https://github.com/gunvirranu/doctest"
    GIT_TAG "v2.4.12-gunvirranu"
)
FetchContent_MakeAvailable(doctest)

add_executable(test_magneto test_magneto.cpp)
target_link_libraries(test_magneto PRIVATE magneto doctest)

# Restore previous
if(DEFINED CMAKE_CXX_CLANG_TIDY_save)
    set(CMAKE_CXX_CLANG_TIDY "${CMAKE_CXX_CLANG_TIDY_save}")
endif()

target_compile_features(test_magneto PRIVATE cxx_std_11)
set_target_properties(
    magneto
    PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED YES
    C_EXTENSIONS NO
)

enable_testing()
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(test_magneto)

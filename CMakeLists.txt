cmake_minimum_required(VERSION 3.14)

# ---- In-source guard ----
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported")
endif()

# ---- Project -----
project(
    magneto
    VERSION 0.0.0
    DESCRIPTION "A portable C library for geomagnetic models"
    HOMEPAGE_URL "https://github.com/gunvirranu/magneto"
    LANGUAGES C
)

option(
    magneto_SINGLE_PRECISION_FLOAT
    "Use single-precision floating point instead of the default double"
    OFF
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    magneto
    src/magneto.c src/model.c
)

target_include_directories(
    magneto ${warning_guard}
    PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(magneto PUBLIC c_std_99)

if(magneto_SINGLE_PRECISION_FLOAT)
    target_compile_definitions(magneto PUBLIC MAGNETO_SINGLE_PRECISION)
endif()

# ---- Developer mode ----
if(NOT magneto_DEVELOPER_MODE)
    return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
    message(AUTHOR_WARNING "Developer mode enabled")
endif()
# Anything below this is purely for developer mode

include(CTest)
if(BUILD_TESTING)
 add_subdirectory(tests)
endif()

cmake_minimum_required(VERSION 3.16)
project(DSP_Guide C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Build flags
add_compile_options(-Wall -Wextra -g -O1)

include_directories(include)

# -- Collect all implementation sources --------------------------------------
file(GLOB_RECURSE DSP_SOURCES
    src/transforms/*.c
    src/filtering/*.c
    src/operations/*.c
    src/spectral/*.c
    src/sampling/*.c
    src/wavelet/*.c
    src/coding/*.c
    src/adaptive/*.c
    src/modulation/*.c
    src/image/*.c
)

# -- Shared library (all DSP implementations) ---------------------------------
add_library(dsp STATIC ${DSP_SOURCES})
target_include_directories(dsp PUBLIC include)
# libm is needed for sin/cos/cexp etc.
target_link_libraries(dsp m)

# -- Demo runner --------------------------------------------------------------
add_executable(dsp_demo src/main.c)
target_link_libraries(dsp_demo dsp m)

# -- Tests --------------------------------------------------------------------
enable_testing()
file(GLOB TEST_SOURCES tests/*.c)
foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} dsp m)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

message(STATUS "DSP Guide -- ${CMAKE_BUILD_TYPE} build")
message(STATUS "Sources found: ${DSP_SOURCES}")

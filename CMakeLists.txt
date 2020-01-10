cmake_minimum_required(VERSION 3.16)

project (snow)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(PROJECT_SOURCE_DIRECTORY ${CMAKE_BINARY_DIR}/src)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# GTest
ADD_SUBDIRECTORY(extlib/gtest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Global include dirs both test and src
include_directories(${CMAKE_BINARY_DIR}/extlib/glm)
include_directories(${PROJECT_SOURCE_DIRECTORY}/solver)
include_directories(${PROJECT_SOURCE_DIRECTORY}/utils)

# Testing
ADD_SUBDIRECTORY(test)

# Binaries
add_subdirectory(src)

cmake_minimum_required(VERSION 3.16)

project (snow)

set(CMAKE_BUILD_TYPE Release)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

message(${MSVC_VERSION})

set(PROJECT_SOURCE_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# GTest
ADD_SUBDIRECTORY(extlib/gtest)
ADD_SUBDIRECTORY(extlib/svd)

# Todo:  The VCPKG doesn't support a cmake FindOpenVdb, and the integrate function doesn't seem to be 
# working with msbuild for me, so just manually specify the paths for now

include_directories(${_VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include)
link_directories(${_VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib)

message(${_VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib)

enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Global include dirs both test and src
include_directories(${CMAKE_SOURCE_DIR}/extlib/glm)
include_directories(${CMAKE_SOURCE_DIR}/extlib/svd)
include_directories(${CMAKE_SOURCE_DIR}/extlib/eigen/Eigen)
include_directories(${CMAKE_SOURCE_DIR}/extlib/openvdb-7.0.0/openvdb)
include_directories(${PROJECT_SOURCE_DIRECTORY}/solver)
include_directories(${PROJECT_SOURCE_DIRECTORY}/utils)
include_directories(${PROJECT_SOURCE_DIRECTORY}/utils)


# Testing
ADD_SUBDIRECTORY(test)

# Binaries
ADD_SUBDIRECTORY(src)

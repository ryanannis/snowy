project(utils)

add_library(utils "")

target_sources(
    utils
    PRIVATE
        Common.hpp
        Math.cpp
        Math.hpp
    PUBLIC
)

target_include_directories(
    solverlib
    PRIVATE
    ../../../extlib/glm
    ../../utils
    PUBLIC
)
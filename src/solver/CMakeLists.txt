project(solver1)

add_library(solverlib "")

target_sources(
    solverlib
    PRIVATE
        CPUSolver.cpp
        CPUSolver.hpp
    PUBLIC
)

target_include_directories(
    solverlib
    PRIVATE
    ../../extlib/glm
    ../utils
    PUBLIC
)
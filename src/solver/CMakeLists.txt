project(solver1)

add_library(solverlib "")

target_sources(
    solverlib
    PRIVATE
        CPUSolver.cpp
        CPUSolver.hpp
        Grid.cpp
        Grid.hpp
        Particle.cpp
        Particle.hpp
    PUBLIC
)

target_include_directories(
    solverlib
    PRIVATE
    ../../extlib/glm
    ../utils
    PUBLIC
)

target_link_libraries(
    solverlib
    utils
)

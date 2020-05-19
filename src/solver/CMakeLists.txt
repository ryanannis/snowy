project(solver1)

add_library(solverlib "")

target_sources(
    solverlib
    PRIVATE
        CPUSolver.cpp
        CPUSolver.hpp
        Grid.cpp
        Grid.hpp
        ParticleSystem.cpp
        ParticleSystem.hpp
        SimulationParameters.cpp
        SimulationParameters.hpp
        SimulationOutput.hpp
        SimulationOutput.cpp
        OvdbConverter.hpp
        OvdbConverter.cpp
        Solver.hpp
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

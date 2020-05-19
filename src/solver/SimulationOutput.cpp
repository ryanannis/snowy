#include "SimulationOutput.hpp"

SimulationOutput::SimulationOutput(const std::vector<Particle>& particles) :
    mParticles(particles)
{
}

const std::vector<Particle>& SimulationOutput::GetParticles() const
{
    return mParticles;
}
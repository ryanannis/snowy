#pragma once

#include "Common.hpp"

// For Particle class
// todo:  Provide a better interface that doesn't expose implementation details
#include "ParticleSystem.hpp"

class SimulationOutput
{
public:
    SimulationOutput(const std::vector<Particle>& particles);
    const std::vector<Particle>& GetParticles() const; // this is nasty shit this class should just have accessors for the particles...

private:
    const std::vector<Particle> mParticles;
};

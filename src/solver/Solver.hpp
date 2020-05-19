#pragma once

#include "Common.hpp"

#include <vector>
#include <memory>

#include "SimulationOutput.hpp"

class Grid;
class ParticleSystem;

class Solver
{
public:
    Solver() = default;

    virtual void AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass) = 0;

    // Each time this is called - the particle list from last frame is invalidated
    virtual void NextFrame() = 0;
    virtual const std::shared_ptr<SimulationOutput> GetOutput() = 0;
};

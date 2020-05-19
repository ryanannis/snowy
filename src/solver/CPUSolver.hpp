#pragma once

#include "Common.hpp"
#include "Solver.hpp"

#include "ParticleSystem.hpp"
#include "Grid.hpp"

class Grid;
class ParticleSystem;

class CPUSolver
{
public:
    CPUSolver(const IVec3& gridDimensions, Float frameLength, const SimulationParameters& params);

    virtual void AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass);

    virtual void NextFrame();

    // Returns the list of particles present in the current simulation step of this solver.
    // Performs a full copy of the particle list so only call this when needed.
    virtual const std::shared_ptr<SimulationOutput> GetOutput();

private:
    void Step(Float timestep);

    Float mFrameLength;
    std::unique_ptr<Grid> mGrid;
    std::unique_ptr<ParticleSystem> mParticleSystem;
    Uint mStepNum;
};
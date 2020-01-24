#pragma once

#include "Common.hpp"

#include <vector>
#include <memory>

#include "ParticleSystem.hpp"
#include "Grid.hpp"


class Grid;
class ParticleSystem;

class CPUSolver
{
public:
    CPUSolver(const IVec3& gridDimensions, Float frameLength, const SimulationParameters& params);

    void AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass);

    // Each time this is called - the particle list from last frame is invalidated
    const std::vector<Particle>& NextFrame();

private:
    void Step(Float timestep);

    Float mFrameLength;
    std::unique_ptr<Grid> mGrid;
    std::unique_ptr<ParticleSystem> mParticleSystem;
    Uint mStepNum;
};



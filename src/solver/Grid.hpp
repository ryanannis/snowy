#pragma once

#include "Common.hpp"

#include "Math.hpp"
#include <iostream>
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include "SimulationParameters.hpp"
#include <mutex>

// todo:  this is needed in here because we have the Weighting templates...  we should just move those somewhere else...
#include "ParticleSystem.hpp"

class Grid;
class MTIterator;
class ParticleSystem;

std::ostream &operator<<(std::ostream &os, Grid const &g);

// Accepts a lambda (IVec grid_coordinate, Float weight)
// BIG TODO:
// This is an optimization nightmare
// (1) Verify if the lambda call optimizes automatically with gdb and MSVC, and if there's any way to do that.
//     This is the hottest function in the program.
// (2) This currently iterates over a radius of 3 even though the kernel radius of 2.  This was because there's
//     some edge cases involving floats equaling integer values - and I didn't want to take any chances before
//     I could get some correct verification data.  I should come back to this later (also possibly cache weight
//     values?)
//
template<typename Func>
void WeightOverParticleNeighbourhood(const SimulationParameters& s, const Particle& p, Func f)
{
    for (Uint i = 0; i < p.numNeighbours; i++)
    {
        f(p.neighbours_coords[i], p.neighbours_nx[i]);
    }
}

template<typename Func>
void WeightGradOverParticleNeighbourhood(const SimulationParameters& s, const Particle& p, Func f)
{
    for (Uint i = 0; i < p.numNeighbours; i++)
    {
        f(p.neighbours_coords[i], p.neighbours_nxgrad[i]);
    }
}

class Cell {
public:
    Cell() = default;

    // We manually lock and unlock mutexes per Cell because 
    // we might not to lock everywhere and operations involving
    // cells are usually very performance sensitive
    inline void Lock() {
        mMutex.lock();
    }
    
    inline void Unlock() 
    {
        mMutex.unlock();
    }

    // Inputs - transferred from particles each step
    Float Mass = 0;
    Vec3 Velocity = {};

    // Intermediate calculation
    Vec3 Force = {};

    // Transferred back to particles each step
    Vec3 VelocityStar = {};
    Vec3 VelocityNext = {};

private:
    std::mutex mMutex;
};

class Grid {
public:
    Grid(const SimulationParameters& params, const IVec3& dims);

    Cell& Get(Uint i, Uint j, Uint k);
    const Cell& Get(Uint i, Uint j, Uint k) const;

    const IVec3& Dims() const;

    void RasterizeParticlesToGrid(const ParticleSystem& ps, MTIterator& mt);
    void ComputeGridForces(const ParticleSystem& ps, MTIterator& mt);
    void UpdateGridVelocities(Float timestep, MTIterator& mt);
    void DoGridBasedCollisions(Float timestep, MTIterator& mt);
    void SolveLinearSystem(Float timestep, MTIterator& mt);
    void ResetGrid();

private:
    Uint coordToIdx(Uint i, Uint j, Uint k) const;

    const SimulationParameters mParams;
    const IVec3 mDims;

    std::vector<Cell> mCells;
    
    friend std::ostream &operator<<(std::ostream &os, Grid const &g);
};

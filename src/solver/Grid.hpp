#pragma once

#include "Common.hpp"

#include "Math.hpp"
#include <iostream>
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include "SimulationParameters.hpp"
#include <mutex>

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
void WeightOverParticleNeighbourhood(const SimulationParameters& s, const Grid& g, const Vec3& Position, Func f) {
    for (int i = -2; i < 3; i++)
    {
        for (int j = -2; j < 3; j++)
        {
            for (int k = -2; k < 3; k++)
            {
                // The coordinate of the cell we are rasterizing to
                int ix = static_cast<int>(Position.x / s.H) + i;
                int iy = static_cast<int>(Position.y / s.H) + j;
                int iz = static_cast<int>(Position.z / s.H) + k;

                // Check bounds
                if (ix < 0 || iy < 0 || iz < 0 || 
                    ix >= g.Dims().x || iy >= g.Dims().y || iz>= g.Dims().z) {
                    continue;
                }

                Float nx = gridWeight(s.H, Position.x / s.H, ix, Position.y / s.H, iy, Position.z / s.H, iz);

                if (nx == 0) {
                    // PERF TODO:  we're going over a 5*5 instead of 4*4 to avoid stupid boundary cases with inequalities
                    // then skipping when the grad is 0...  we can do better than this
                    continue;
                }

                f(IVec3(ix, iy, iz), nx);
            }
        }
    }
}

template<typename Func>
void WeightGradOverParticleNeighbourhood(const SimulationParameters& s, const Grid& g, const Vec3& Position, Func f) {
    for (int i = -2; i < 3; i++)
    {
        for (int j = -2; j < 3; j++)
        {
            for (int k = -2; k < 3; k++)
            {
                // The coordinate of the cell we are rasterizing to
                int ix = static_cast<int>(Position.x / s.H) + i;
                int iy = static_cast<int>(Position.y / s.H) + j;
                int iz = static_cast<int>(Position.z / s.H) + k;

                // Check bounds
                if (ix < 0 || iy < 0 || iz < 0 ||
                    ix >= g.Dims().x || iy >= g.Dims().y || iz >= g.Dims().z) {
                    continue;
                }

                Vec3 nxgrad = gridWeightGrad(s.H, Position.x / s.H, ix, Position.y / s.H, iy, Position.z / s.H, iz);
                Float nx = gridWeight(s.H, Position.x / s.H, ix, Position.y / s.H, iy, Position.z / s.H, iz);
                if (nxgrad != Vec3(0))
                {
                    continue;
                    // PERF TODO:  we're going over a 5*5 instead of 4*4 to avoid stupid boundary cases with inequalities
                    // then skipping when the grad is 0...  we can do better than this
                }

                f(IVec3(ix, iy, iz), nxgrad);
            }
        }
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

#pragma once

#include "Common.hpp"

#include <iostream>
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include "SimulationParameters.hpp"

class Grid;
class ParticleSystem;

std::ostream &operator<<(std::ostream &os, Grid const &g);

// Yeah this shouldn't be here
Vec3 gridWeightGrad(const SimulationParameters& s, Float x, Float ix, Float y, Float iy, Float z, Float iz);
Float gridWeight(const SimulationParameters& s, Float x, Float ix, Float y, Float iy, Float z, Float iz);

// Accepts a lambda (IVec grid_coordinate, Float weight)
template<typename Func>
void WeightOverParticleNeighbourhood(const SimulationParameters& s, const Grid& g, const Vec3& Position, Func f) {
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            for (int k = -1; k < 2; k++)
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

                Float nx = gridWeight(s, Position.x / s.H, ix, Position.y / s.H, iy, Position.z / s.H, iz);

                f(IVec3(ix, iy, iz), nx);
            }
        }
    }
}

template<typename Func>
void WeightGradOverParticleNeighbourhood(const SimulationParameters& s, const Grid& g, const Vec3& Position, Func f) {
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            for (int k = -1; k < 2; k++)
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

                Vec3 nxgrad = gridWeightGrad(s, Position.x / s.H, ix, Position.y / s.H, iy, Position.z / s.H, iz);

                f(IVec3(ix, iy, iz), nxgrad);
            }
        }
    }
}

class Cell {
public:
    Cell() = default;

    // Inputs - transferred from particles each step
    Float Mass = 0;
    Vec3 Velocity = {};

    // Intermediate calculation
    Vec3 Force = {};

    // Transferred back to particles each step
    Vec3 VelocityStar = {};
    Vec3 VelocityNext = {};

};

class Grid {
public:
    Grid(const SimulationParameters& params, const IVec3& dims);

    Cell& Get(Uint i, Uint j, Uint k);
    Cell Get(Uint i, Uint j, Uint k) const;

    const IVec3& Dims() const;

    void RasterizeParticlesToGrid(const ParticleSystem& ps);
    void ComputeGridForces(const ParticleSystem& ps);
    void UpdateGridVelocities(Float timestep);
    void DoGridBasedCollisions(Float timestep);
    void SolveLinearSystem(Float timestep);

private:
    Uint coordToIdx(Uint i, Uint j, Uint k) const;

    const SimulationParameters mParams;
    const IVec3 mDims;

    std::vector<Cell> mCells;
    
    friend std::ostream &operator<<(std::ostream &os, Grid const &g);
};

#pragma once

#include "Common.hpp"

#include <iostream>
#include <vector>
#include <memory>
#include <glm/glm.hpp>

class Grid;
class Particle;

std::ostream &operator<<(std::ostream &os, Grid const &g);

// Yeah this shouldn't be here
Vec3 gridWeightGrad(Float x, Float ix, Float y, Float iy, Float z, Float iz);
Float gridWeight(Float x, Float ix, Float y, Float iy, Float z, Float iz);

class Cell {
public:
    Cell() = default;

    Float Mass = 0;
    Vec3 Velocity = {};
    Vec3 VelocityStar = {};
    Vec3 Force = {};
};

class Grid {
public:
    Grid(Uint lenI, Uint lenJ, Uint leK);

    // Returns the cell at the given instance
    // For reducing bounds checking elsewhere - will return an empty cell if accessed at an invalid index

    // Todo:  It has a performance overhead to do the bounds checking here
    // Todo:  Given my operations - this is likely not the most cache-efficent way to split up the cell parameters
    Cell& Get(Uint i, Uint j, Uint k);
    Cell Get(Uint i, Uint j, Uint k) const;

    Uint LenI() const;
    Uint LenJ() const;
    Uint LenK() const;

    void TransferParticlesToGrid(std::vector<Particle>& particles);
    void EstimateParticleVolumes(std::vector<Particle>& particles);
    void ComputeGridForces(std::vector<Particle>& particles);
    void UpdateGridVelocities(Float timestep);
    void DoGridBasedCollisions(Float timestep);
    void VelocityImplicitUpdate(Float timestep);

private:
    std::vector<std::vector<Uint>> CreateCellToParticleMap(std::vector<Particle>& particles) const;

    const Uint mLenI;
    const Uint mLenJ;
    const Uint mLenK;

    Uint coordToIdx(Uint i, Uint j, Uint k) const;

    // Used for the mechanism of returning an "empty cell"
    Cell mDummyCell;
    std::vector<Cell> mCells;
    
    friend std::ostream &operator<<(std::ostream &os, Grid const &g);
};

#pragma once

#include "Common.hpp"

#include <iostream>
#include <vector>
#include <memory>
#include <glm/glm.hpp>

class Particle {
public:
    Particle(const glm::vec3& pos);
    Float& Mass();
    Float Mass() const;
    const glm::vec3 Position() const;

private:
    glm::vec3 mPos;
    Float mMass = 0;
};

class Cell {
public:
    Cell() = default;
    Float& Mass();
    Float Mass() const;


private:
    Float mMass = 0;
};


class Grid;

std::ostream &operator<<(std::ostream &os, Grid const &g);

class Grid {
public:
    Grid(Uint lenI, Uint lenJ, Uint leK);

    // Returns the cell at the given instance
    // For reducing bounds checking elsewhere - will return an empty cell if accessed at an invalid index

    // Todo:  It has a performance overhead to do the bounds checking here
    // Todo:  Given my operations - this is likely not the most cache-efficent way to split up the cell parameters
    Cell& Get(Uint i, Uint j, Uint k);
    Cell Get(Uint i, Uint j, Uint k) const;

private:
    const Uint mLenI;
    const Uint mLenJ;
    const Uint mLenK;

    // Used for the mechanism of returning an "empty cell"
    Cell mDummyCell;
    std::vector<Cell> mCells;
    
    friend std::ostream &operator<<(std::ostream &os, Grid const &g);
};

class CPUSolver
{
public:
    CPUSolver();
    void Step(Uint substeps);

private:

    void RasterizeParticlesToGrid();
    void ComputeParticleVolumesAndDensities();
    
    std::unique_ptr<Grid> mGrid;
    std::vector<Particle> mParticles;
};


// todo:  this is uh questionable design
class RasterizationUtils
{
public:
    static void RasterizeMassToGrid(Grid& g, const std::vector<Particle>& particles);
};
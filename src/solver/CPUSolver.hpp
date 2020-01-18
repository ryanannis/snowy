#pragma once

#include "Common.hpp"

#include <memory>

#include "Particle.hpp"
#include "Grid.hpp"


class CPUSolver
{
public:
    CPUSolver();
    void Step(Uint substeps);

private:
    std::unique_ptr<Grid> mGrid;
    std::vector<Particle> mParticles;
};



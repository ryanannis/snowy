#include "Grid.hpp"
#include "ParticleSystem.hpp"
#include "Multithread.hpp"

// This is called with coordinates normalized to the size of the grid
// (eg. divided by h)

Grid::Grid(const SimulationParameters& params, const IVec3& dims) :
    mParams(params),
    mDims(dims),
    mCells(dims.x * dims.y * dims.z)
{
} 

Cell& Grid::Get(Uint i, Uint j, Uint k)
{
    return mCells[coordToIdx(i, j, k)];
}

const Cell& Grid::Get(Uint i, Uint j, Uint k) const
{
    return mCells[coordToIdx(i, j, k)];
}

const IVec3& Grid::Dims() const
{
    return mDims;
}

Uint Grid::coordToIdx(Uint i, Uint j, Uint k) const
{
    return i + mDims.x * j + mDims.x * mDims.y * k;
}

void Grid::RasterizeParticlesToGrid(const ParticleSystem& ps, MTIterator& mt)
{
    mt.IterateOverVector(ps.GetParticles(), [&](const Particle& particle) {
        WeightOverParticleNeighbourhood(mParams, particle,
            [&](IVec3 pos, Float weight) {
                // Transfer mass
                Cell& c = Get(pos.x, pos.y, pos.z);

                c.Lock();
                c.Mass += weight * particle.mass;

                // Transfer velocity (normalized)
                c.Velocity += particle.velocity * particle.mass * weight;
                c.Unlock();
            });
    });
}

void Grid::ComputeGridForces(const ParticleSystem& ps, MTIterator& mt)
{
    mt.IterateOverVector(ps.GetParticles(), [&](const Particle& particle) {
        WeightGradOverParticleNeighbourhood(mParams, particle,
            [&](IVec3 pos, Vec3 weightgrad) {
                Mat3 stress = -ps.CalculateCauchyStress(particle);
                Vec3 dforce = particle.volume * stress * weightgrad;
                ASSERT_VALID_VEC3(dforce);
                
                Cell& c = Get(pos.x, pos.y, pos.z);
                c.Lock();
                c.Force += dforce;
                c.Unlock();
            }
        );
    });
}

void Grid::UpdateGridVelocities(Float timestep, MTIterator& mt) {
    mt.IterateOverVector(mCells, [&](Cell& c) {
        if(c.Mass > 0) {
            c.Velocity /= c.Mass; // normalize velocity for energy conservation
            c.Force += Vec3(0.0, 0.0, mParams.GRAVITY * c.Mass);
            c.VelocityStar += c.Velocity + timestep / c.Mass * c.Force;
            ASSERT_VALID_VEC3(c.VelocityStar);
        }
    });
}

void Grid::DoGridBasedCollisions(Float timestep, MTIterator& mt)
{
    // todo:
}

void Grid::SolveLinearSystem(Float timestep, MTIterator& mt)
{
    mt.IterateOverVector(mCells, [&](Cell& c) {
        c.VelocityNext = c.VelocityStar;   
    });
}

void Grid::ResetGrid()
{
    for (Cell& c : mCells)
    {
        c.Force = Vec3(0.0);
        c.Mass = 0.0;
        c.Velocity = Vec3(0.0);
        c.VelocityNext = Vec3(0.0);
        c.VelocityStar = Vec3(0.0);
    }
}

std::ostream &operator<<(std::ostream &os, Grid const &g) { 
    for (int i = 0; i < g.mDims.x; i++)
    {
        for (int j = 0; j < g.mDims.y; j++)
        {
            for (int k = 0; k < g.mDims.z; k++)
            {
                os << g.Get(i, j, k).Mass << " ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    return os;
}

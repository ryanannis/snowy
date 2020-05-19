#include "Grid.hpp"
#include "ParticleSystem.hpp"

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

Cell Grid::Get(Uint i, Uint j, Uint k) const
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

void Grid::RasterizeParticlesToGrid(const ParticleSystem& ps)
{
    for (const auto& particle : ps.GetParticles())
    {
        WeightOverParticleNeighbourhood(
            mParams, *this, particle.pos,
            [&](IVec3 pos, Float weight) {
                // Transfer mass
                Cell& c = Get(pos.x, pos.y, pos.z);
                Get(pos.x, pos.y, pos.z).Mass += weight * particle.mass;

                // Transfer velocity (normalized)
                Get(pos.x, pos.y, pos.z).Velocity += particle.velocity * particle.mass * weight;
            }
        );
    }
}

void Grid::ComputeGridForces(const ParticleSystem& ps)
{
    for (const auto& particle : ps.GetParticles()) {
        WeightGradOverParticleNeighbourhood(
            mParams, *this, particle.pos,
            [&](IVec3 pos, Vec3 weightgrad) {
                Mat3 stress = -ps.CalculateCauchyStress(particle);
                Vec3 dforce = particle.volume * stress * weightgrad;
                ASSERT_VALID_VEC3(dforce);

                Get(pos.x, pos.y, pos.z).Force += dforce;
                
            }
        );
    }
}

void Grid::UpdateGridVelocities(Float timestep) {
    for (auto& c : mCells) {        
        if(c.Mass > 0) {
            c.Velocity /= c.Mass; // normalize velocity for energy conservation
            c.Force += Vec3(0.0, 0.0, mParams.GRAVITY * c.Mass);
            c.VelocityStar += c.Velocity + timestep / c.Mass * c.Force;
            ASSERT_VALID_VEC3(c.VelocityStar);
        }
    }
}

void Grid::DoGridBasedCollisions(Float timestep)
{
    // todo:
}

void Grid::SolveLinearSystem(Float timestep)
{
    // todo:  only explicit update right nowD
    for (auto& c : mCells) {
        c.VelocityNext = c.VelocityStar;
    }
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

#include "Grid.hpp"
#include "ParticleSystem.hpp"

Float N_x(Float x)
{
    Float res;

    Float x2 = x * x;
    Float abx = std::abs(x);
    Float x3 = abx * x2;

    // Check we didn't feed any cells outside kernel boundaries
    assert(abx <= 2.0);

    if (abx < 1)
    {
        res =  (0.5 * x3) - (x2) + (2.0 / 3.0); 
    }
    else
    {
        res = (-1.0 / 6.0 * x3) + (x2) - (2.0 * abx) + (4.0/3.0); 
    }

    ASSERT_VALID_FLOAT(res);

    return res;
}


Float dN_x(Float x)
{
    Float res;

    Float x2 = x * x;
    Float abx = std::abs(x);
    Float sgnx = x / std::abs(x);

    // Check we didn't feed any cells outside kernel boundaries
    assert(abx <= 2.0);

    if (abx < EPSILON) {
        res = 0;
    }
    else if (abx < 1)
    {
        res =  (1.5 * x2 * sgnx) - (2.0 * x); 
    }
    else
    {
        res = (-0.5 * x2 * sgnx) + (2.0 * x) - (2.0 * sgnx); 
    }

    ASSERT_VALID_FLOAT(res);

    return res;
}

// This is called with coordinates normalized to the size of the grid
// (eg. divided by h)
Float gridWeight(const SimulationParameters& s, Float x, Float ix, Float y, Float iy, Float z, Float iz)
{
    return N_x(x - ix) * N_x(y - iy) * N_x(z - iz);
}

Vec3 gridWeightGrad(const SimulationParameters& s, Float x, Float ix, Float y, Float iy, Float z, Float iz)
{
    Float dx = x - ix;
    Float dy = y - iy;
    Float dz = z - iz;

    return Vec3(
        dN_x(dx) * N_x(dy)  * N_x(dz),
        N_x(dx)  * dN_x(dy) * N_x(dz),
        N_x(dx)  * N_x(dy)  * dN_x(dz)
    ) / s.H; // todo:  is this derivative right??
}


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
    return i * mDims.x * mDims.y + j * mDims.y + k;
}

void Grid::RasterizeParticlesToGrid(const ParticleSystem& ps)
{
    for (const auto& particle : ps.GetParticles())
    {
        WeightOverParticleNeighbourhood(
            mParams, *this, particle.pos,
            [&](IVec3 pos, Float weight) {
                // Transfer mass
                Get(pos.x, pos.y, pos.z).Mass += weight * particle.mass;

                // Transfer velocity (normalized)
                Get(pos.x, pos.y, pos.z).Velocity += particle.velocity * particle.mass * weight;
            }
        );
    }

    // Normalize
    for (auto& cell : mCells) {
        if (cell.Mass == 0) continue;
        cell.Velocity /= cell.Mass;
    }
}

void Grid::ComputeGridForces(const ParticleSystem& ps)
{
    for (const auto& particle : ps.GetParticles()) {
        WeightGradOverParticleNeighbourhood(
            mParams, *this, particle.pos,
            [&](IVec3 pos, Vec3 weightgrad) {
                Mat3 stress = ps.CalculateCauchyStress(particle);
                Vec3 dforce = particle.volume * stress * weightgrad;
                ASSERT_VALID_VEC3(dforce);

                Get(pos.x, pos.y, pos.z).Force += dforce;
            }
        );
    }
}

void Grid::UpdateGridVelocities(Float timestep) {
    for (auto& c : mCells) {
        // Add gravity
        c.Force += Vec3(0.0, 0.0, mParams.GRAVITY * c.Mass);

        c.Velocity = c.VelocityStar;
        
        // Possible to be zero if a particle moved out of the cell in the last time step
        if(c.Mass > 0) {
            c.VelocityStar += timestep / c.Mass * c.Force;
            ASSERT_VALID_VEC3(c.VelocityStar);
        }
        
        // Reset
        c.Force = Vec3(0.0);
        c.Mass = 0;
        c.Velocity = Vec3(0.0);
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

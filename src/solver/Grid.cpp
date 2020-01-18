#include "Grid.hpp"
#include "Particle.hpp"

Float N_x(Float x)
{
    Float res;

    Float x2 = x * x;
    Float abx = std::abs(x);
    Float x3 = abx * x2;

    assert(abx < 2.5);

    if (abx < 1)
    {
        res =  (0.5 * x3) - (x2) + (2.0 / 3.0); 
    }
    else
    {
        res = (-1.0 / 6.0 * x3) + (x2) - (2.0 * abx) + (4.0/3.0); 
    }

    assert(res > 0.0);
    assert(res < 1.0);

    return res / H;
}


Float dN_x(Float x)
{
    Float res;

    Float x2 = x * x;
    Float abx = std::abs(x);
    Float sgnx = x / std::abs(x);

    assert(abx < 2.5);

    if (abx < 1)
    {
        res =  (1.5 * x2 * sgnx) - (2 * x); 
    }
    else
    {
        res = (-0.5 * x2 * sgnx) + (2.0 * x) - (2.0 * sgnx); 
    }

    assert(res > 0.0);
    assert(res < 1.0);

    return res / H;
}

Float gridWeight(Float x, Float ix, Float y, Float iy, Float z, Float iz)
{
    return N_x(x - ix) * N_x(y - iy) + N_x(z - iz);
}

Vec3 gridWeightGrad(Float x, Float ix, Float y, Float iy, Float z, Float iz)
{
    Float dx = x - ix;
    Float dy = y - iy;
    Float dz = z - iz;

    return Vec3(
        dN_x(dx) * N_x(dy) + N_x(dz),
        N_x(dx) * dN_x(dy) + N_x(dz),
        N_x(dx) * N_x(dy) + dN_x(dz)
    ); // todo:  the H?
}


Grid::Grid(Uint lenI, Uint lenJ, Uint lenK) : 
    mLenI(lenI),
    mLenJ(lenJ),
    mLenK(lenK),
    mCells(lenI * lenJ * lenK)
{
}

Cell& Grid::Get(Uint i, Uint j, Uint k)
{
    if (i < 0 || i >= mLenI || 
        j < 0 || j >= mLenJ ||
        k < 0 || k >= mLenK)
    {
        mDummyCell = Cell(); // reset dummycell
        return mDummyCell;
    }

    return mCells[coordToIdx(i, j, k)];
}

Cell Grid::Get(Uint i, Uint j, Uint k) const
{
    if (i < 0 || i >= mLenI || 
        j < 0 || j >= mLenJ ||
        k < 0 || k >= mLenK)
    {
        return Cell();
    }

    return mCells[coordToIdx(i, j, k)];
}

Uint Grid::coordToIdx(Uint i, Uint j, Uint k) const
{
    return i * mLenI * mLenJ + j * mLenJ + k;
}

void Grid::TransferParticlesToGrid(std::vector<Particle>& particles)
{
    // Calculate interpolation weights for a 2x2 kernel
    // todo:  this is repeated many time, but I also don't trust the 
    // compiler enough to use a lambda...
    for (auto& particle : particles) 
    {
        particle.Mass = 0;
        particle.Velocity = {};

        Float x = particle.Position.x;
        Float y = particle.Position.y;
        Float z = particle.Position.z;
        for (int i = -1; i <= 2; i++)
        {
            for (int j = -1; j <= 2; j++)
            {
                for (int k = -1; k <= 2; k++)
                {
                    // The coordinate of the cell we are rasterizing to
                    int ix = static_cast<int>(x) + i;
                    int iy = static_cast<int>(y) + j;
                    int iz = static_cast<int>(z) + k ;

                    Float nx = gridWeight(x, ix, y, iy, z, iz);

                    // Transfer mass
                    Get(ix, iy, iz).Mass += nx * particle.Mass;

                    // Transfer velocity (normalized)
                    Get(ix, iy, iz).Velocity += particle.Velocity * particle.Mass * nx;
                }
            }
        }
    }

    // Iterate over all cells for velocity normalization step
    for (int i = 0; i < LenI(); i++)
    {
        for (int j = 0; j < LenJ(); j++)
        {
            for (int k = 0; k < LenK(); k++)
            {
                Get(i, j, k).Velocity /= Get(i, j, k).Mass;
            }
        }
    }
}

void Grid::EstimateParticleVolumes(std::vector<Particle>& particles)
{
    for (auto& particle : particles) 
    {
        Float x = particle.Position.x;
        Float y = particle.Position.y;
        Float z = particle.Position.z;

        Float particleDensity = 0;

        for (int i = -1; i <= 2; i++)
        {
            for (int j = -1; j <= 2; j++)
            {
                for (int k = -1; k <= 2; k++)
                {
                    // The coordinate of the cell we are rasterizing to
                    int ix = static_cast<int>(x) + i;
                    int iy = static_cast<int>(y) + j;
                    int iz = static_cast<int>(z) + k ;

                    Float cellvolume = H * H * H;

                    Float nx = gridWeight(x, ix, y, iy, z, iz);

                    particleDensity += nx * Get(ix, iy, iz).Mass / cellvolume;
                }
            }
        }

        particle.Volume = particle.Mass / particleDensity ;
    }
}

void Grid::ComputeGridForces(std::vector<Particle>& particles)
{
    // Clear forces
    for (auto& c : mCells) {
        c.Force = {};
    }

    for (auto& particle : particles) 
    {
        Float x = particle.Position.x;
        Float y = particle.Position.y;
        Float z = particle.Position.z;

        for (int i = -1; i <= 2; i++)
        {
            for (int j = -1; j <= 2; j++)
            {
                for (int k = -1; k <= 2; k++)
                {
                    // The coordinate of the cell we are rasterizing to
                    int ix = static_cast<int>(particle.Position.x) + i;
                    int iy = static_cast<int>(particle.Position.y) + j;
                    int iz = static_cast<int>(particle.Position.z) + k ;

                    Vec3 nxGrad = gridWeightGrad(x, ix, y, iy, z, iz);

                    Get(ix, iy, iz).Force += particle.Volume * particle.CalculateCauchyStress() * nxGrad;
                }
            }
        }
    }
}

void Grid::UpdateGridVelocities(Float timestep) {
    for (auto& c : mCells) {
        c.Velocity = c.VelocityStar;
        c.VelocityStar += timestep / c.Mass * c.Force;
    }
}

void Grid::DoGridBasedCollisions(Float timestep)
{
    // todo:
}

void Grid::VelocityImplicitUpdate(Float timestep)
{
    // todo:  only explicit update right nowD
}

std::ostream &operator<<(std::ostream &os, Grid const &g) { 
    for (int i = 0; i < g.mLenI; i++)
    {
        for (int j = 0; j < g.mLenJ; j++)
        {
            for (int k = 0; k < g.mLenK; k++)
            {
                os << g.Get(i, j, k).Mass << " ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    return os;
}

#include "CPUSolver.hpp"
#include "Math.hpp"

#include <assert.h> 
#include <cmath> 

Float CubicRasterizationWeight(Float x)
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

    std::cout << res << std::endl;
    std::cout << x << std::endl;

    assert(res > 0.0);
    assert(res < 1.0);

    return res;
}


Particle::Particle(const glm::vec3& pos) :
    mPos(pos)
{
}

Float& Particle::Mass()
{
    return mMass;
}

Float Particle::Mass() const
{
    return mMass;
}

const glm::vec3 Particle::Position() const
{
    return mPos;
}

Float& Cell::Mass()
{
    return mMass;
}

Float Cell::Mass() const
{
    return mMass;
}

void RasterizationUtils::RasterizeMassToGrid(Grid& g, const std::vector<Particle>& particles)
{
    // Calculate interpolation weights for a 2x2 kernel
    // todo:  this is coded very stupidly - this method is going to be a high
    // percentage of compute time
    for (auto const& particle : particles) 
    {
        Float x = particle.Position().x;
        Float y = particle.Position().y;
        Float z = particle.Position().z;
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

                    Float nx = CubicRasterizationWeight(x - static_cast<Float>(ix)) * 
                                CubicRasterizationWeight(y - static_cast<Float>(iy)) * 
                                CubicRasterizationWeight(z - static_cast<Float>(iz));

                    g.Get(ix, iy, iz).Mass() += nx * particle.Mass();
                }
            }
        }
    }
}

std::ostream &operator<<(std::ostream &os, Grid const &g) { 
    for (int i = 0; i < g.mLenI; i++)
    {
        for (int j = 0; j < g.mLenJ; j++)
        {
            for (int k = 0; k < g.mLenK; k++)
            {
                os << g.Get(i, j, k).Mass() << " ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    return os;
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

    return mCells[i * mLenI * mLenJ + j * mLenJ + k];
}

Cell Grid::Get(Uint i, Uint j, Uint k) const
{
    if (i < 0 || i >= mLenI || 
        j < 0 || j >= mLenJ ||
        k < 0 || k >= mLenK)
    {
        return Cell();
    }

    return mCells[i * mLenI * mLenJ + j * mLenJ + k];
}

CPUSolver::CPUSolver()
{
}

void CPUSolver::Step(Uint substeps)
{
}

void CPUSolver::RasterizeParticlesToGrid()
{
}

void CPUSolver::ComputeParticleVolumesAndDensities()
{
}
#include "CPUSolver.hpp"
#include "Math.hpp"

#include <assert.h> 
#include <cmath> 

Particle::Particle() 
{
}

Cell::Cell()
{
}

Grid::Grid(Uint lenI, Uint lenJ, Uint lenK) : 
    mLenI(lenI),
    mLenJ(lenJ),
    mLenK(lenK),
    mCells(lenI * lenJ * lenK)
{
}

Cell& Grid::get(Uint i, Uint j, Uint k)
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

CPUSolver::CPUSolver()
{

}

void CPUSolver::Step(Uint substeps)
{

}

double CubicRasterizationWeight(double x)
{
    double res;

    double x2 = x;
    double abx = std::abs(x);
    double x3 = abx * x2;

    if (abx < 1)
    {
        res =  (0.5 * x3) - (x2) + (2.0 / 3.0); 
    }
    else
    {
        res = (-1.0 / 6.0 * x3) + (x2) - (2.0 * abx) + (4.0/3.0); 
    }

    assert(res > 0.0);
    assert(res <= 1.0);

    return res;
}

void CPUSolver::RasterizeParticlesToGrid()
{
    // Calculate interpolation weights for a 2x2 kernel
    // todo:  this is coded very stupidly - this method is going to be a high
    // percentage of compute time
    double val;
    double x;
    double y;
    double z;
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

                double nx = CubicRasterizationWeight(x - static_cast<double>(ix)) * 
                            CubicRasterizationWeight(y - static_cast<double>(iy)) * 
                            CubicRasterizationWeight(z - static_cast<double>(iz));
                

            }
        }
    }
}

void CPUSolver::ComputeParticleVolumesAndDensities()
{

}
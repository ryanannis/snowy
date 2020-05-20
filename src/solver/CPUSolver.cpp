#include "CPUSolver.hpp"

#include <assert.h> 
#include <cmath> 

#include "Math.hpp"


CPUSolver::CPUSolver(const IVec3& gridDimensions, Float frameLength, const SimulationParameters& params) :
    mGrid(std::make_unique<Grid>(params, gridDimensions)),
    mParticleSystem(std::make_unique<ParticleSystem>(params)),
    mFrameLength(frameLength),
    mStepNum(0),
    mMt(12) // TODO:  make customization - 12 threads for now
{
}

void CPUSolver::Step(Float timestep)
{
    std::cout << "Beginning step " << mStepNum << "." << std::endl;

    mParticleSystem->CacheParticleGrads(*mGrid, mMt);

    // @1:  Rasterize particle data to the grid
    mGrid->RasterizeParticlesToGrid(*mParticleSystem, mMt);

    // @2:  Compute particle volumes and densities
    if(mStepNum == 0)
        mParticleSystem->EstimateParticleVolumes(*mGrid, mMt);

    // @3: Compute grid forces
    mGrid->ComputeGridForces(*mParticleSystem, mMt);

    // @4: Compute grid forces
    mGrid->UpdateGridVelocities(timestep, mMt);

    // @5:  Grid based body collisions
    mGrid->DoGridBasedCollisions(timestep, mMt);

    // @6:  Solve linear system
    mGrid->SolveLinearSystem(timestep, mMt);

    // @7: Update deformation gradient
    mParticleSystem->UpdateDeformationGradients(timestep, *mGrid, mMt);

    // @8: Update Particle Velocities
    mParticleSystem->UpdateVelocities(*mGrid, mMt);

    // @9: Particle-based body collisions  
    mParticleSystem->BodyCollisions(timestep, mMt);

    // @10:  Update particle positions
    mParticleSystem->UpdatePositions(timestep, mMt);

    // We've transferred everything to the particles.  Reset the accumulators.
    mGrid->ResetGrid();

    mStepNum++;
}

void CPUSolver::NextFrame()
{
    // Just do a constant timestep for now...  (TODO)
    const Float timestep = 0.0001;

    // Framelength needs to be divisible by timestep
    const Uint stepsPerFrame = std::floor(mFrameLength / timestep);
    assert(stepsPerFrame != 0);

    for (Uint i = 0; i < stepsPerFrame; i++) {
        std::cout << "Time: " << Float(mStepNum) * mFrameLength / Float(stepsPerFrame) << std::endl;
        Particle& p = mParticleSystem->GetParticles()[0];
        std::cout << "[" << p.velocity.x << "," << p.velocity.y << "," << p.velocity.z << "]" << std::endl;
        std::cout << "[" << p.pos.x << "," << p.pos.y << "," << p.pos.z << "]" << std::endl;

        std::cout << "[" << p.m_F_p[0][0] << "," << p.m_F_p[1][1] << "," << p.m_F_p[2][2] << "]" << std::endl;

        Step(mFrameLength / Float(stepsPerFrame));
    }
}

const std::shared_ptr<SimulationOutput> CPUSolver::GetOutput()
{
    return std::shared_ptr<SimulationOutput>(
        new SimulationOutput(mParticleSystem->GetParticles())
    );
}

void CPUSolver::AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass)
{
    mParticleSystem->AddParticle(pos, velocity, mass);
}

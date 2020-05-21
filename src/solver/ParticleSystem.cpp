#include "ParticleSystem.hpp"

#include "Grid.hpp"
#include "Math.hpp"
#include "Multithread.hpp"
#include "glm/gtx/matrix_operation.hpp"

#include <exception>

Particle::Particle(const Vec3& pos, Float mass, const Vec3& velocity) :
    pos(pos),
    mass(mass),
    velocity(velocity),
    volume(0.0), // This is set later
    m_F_p(Mat3(1.0)),
    m_F_e(Mat3(1.0)),
    m_R_e(Mat3(1.0))
{
}

ParticleSystem::ParticleSystem(const SimulationParameters& parameters) :
    mParams(parameters)
{
}

void ParticleSystem::AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass)
{
    mParticles.push_back(Particle(pos, mass, velocity)); // Volume is 0 by default
}

Mat3 ParticleSystem::CalculateVelocityGradient(const Particle& p, const Grid& g) const 
{
    Mat3 velGrad = Mat3(Float(0.0));

    WeightGradOverParticleNeighbourhood(
        mParams, p,
        [&](IVec3 pos, Vec3 weightgrad) {
            velGrad += glm::outerProduct(g.Get(pos.x, pos.y, pos.z).VelocityStar, weightgrad);
        }
    );

    return velGrad;
}

Mat3 ParticleSystem::CalculateCauchyStress(const Particle& p) const
{
    Float j_p = glm::determinant(p.m_F_p);
    Float j_e = glm::determinant(p.m_F_e);

    Float hardening = exp(mParams.HARDENING * (1 - j_p));
    Float mu = mParams.MU_0 * hardening;
    Float lambda = mParams.LAMBDA_0 * hardening;

    auto result = Float(2.0) * mu * (p.m_F_e - p.m_R_e) * glm::transpose(p.m_F_e) + Mat3(lambda * (j_e - 1) * j_e);
    
    ASSERT_VALID_MAT3(result);

    return result;
}
 
void ParticleSystem::BodyCollisions(Float dt, MTIterator& mt)
{
}

void ParticleSystem::CacheParticleGrads(const Grid& g, MTIterator& mt)
{
    mt.IterateOverVector(mParticles, [&](Particle& p) {
        const auto& Position = p.pos;
        const auto& H = mParams.H;

        p.numNeighbours = 0;
        
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                for (int k = -2; k < 3; k++)
                {
                    // The coordinate of the cell we are rasterizing to
                    int ix = static_cast<int>(Position.x / H) + i;
                    int iy = static_cast<int>(Position.y / H) + j;
                    int iz = static_cast<int>(Position.z / H) + k;

                    // Check bounds
                    if (ix < 0 || iy < 0 || iz < 0 ||
                        ix >= g.Dims().x || iy >= g.Dims().y || iz >= g.Dims().z) {
                        continue;
                    }

                    Vec3 nxgrad = gridWeightGrad(H, Position.x / H, ix, Position.y / H, iy, Position.z / H, iz);
                    Float nx = gridWeight(H, Position.x / H, ix, Position.y / H, iy, Position.z / H, iz);

                    if (nxgrad == Vec3(0) && nx == 0)
                    {
                        continue;
                    }

                    
                    if(p.numNeighbours > 63) 
                    {
                        throw std::exception("Error in creating gradient kernel cache, crashing due to potential data corruption!");
                    }

                    p.neighbours_coords[p.numNeighbours] = IVec3(ix, iy, iz);
                    p.neighbours_nx[p.numNeighbours] = nx;
                    p.neighbours_nxgrad[p.numNeighbours] = nxgrad;
                    p.numNeighbours++;

                }
            }
        }
    });
}

void ParticleSystem::UpdateDeformationGradients(Float dt, const Grid& g, MTIterator& mt)
{
    mt.IterateOverVector(mParticles, [&](Particle& p) {
        // First attribute all new changes to elastic part of deformation
        Mat3 new_fe = (Mat3(Float(1.0)) + dt * CalculateVelocityGradient(p, g)) * p.m_F_e;
        Mat3 new_f = new_fe * p.m_F_p;
        
        Mat3 u(1.0);
        Mat3 s(1.0);
        Mat3 v(1.0);
        svd3(new_fe, u, s, v);

        // Todo:  Add an assert for reconstructing the SVD - I have a bad history with SVD libraries

        Mat3 sinv(1.0);

        // Clamp the singular values
        for(Uint i = 0; i < 3; i++) {
            s[i][i] = glm::clamp(s[i][i], Float(1.0) - mParams.PHI_C, Float(1.0) + mParams.PHI_S);
            sinv[i][i] = Float(1.0) / s[i][i];
        }

        p.m_F_e = u * s * glm::transpose(v);
        p.m_F_p = v * sinv * glm::transpose(u) * new_f;

        // todo:  can we use any properties of the previous SVD to speed this up?

        svd3(p.m_F_e, u, s, v); 
        p.m_R_e = u * glm::transpose(v);

        ASSERT_VALID_MAT3(p.m_F_e);
        ASSERT_VALID_MAT3(p.m_F_p);
        ASSERT_VALID_MAT3(p.m_R_e);
    });
}

void ParticleSystem::EstimateParticleVolumes(const Grid& g, MTIterator& mt)
{
    mt.IterateOverVector(mParticles, [&](Particle& particle) {
        Float particleDensity = 0;

        WeightOverParticleNeighbourhood(
            mParams, particle,
            [&](IVec3 pos, Float weight) {

                Float cellvolume = mParams.H * mParams.H * mParams.H;
                particleDensity += weight * g.Get(pos.x, pos.y, pos.z).Mass / cellvolume;
            }
        );

        // This ((should)) never be < 0 because the mass is rasterized to the same grid cells
        // the step before this
        ASSERT_VALID_FLOAT(particleDensity);
        if(particleDensity > 0) {
            particle.volume = particle.mass / particleDensity;
        }
    });
}


void ParticleSystem::CalculateFlipPicVelocity(const Particle& p, const Grid& g, Vec3& flip, Vec3& pic) const
{
    flip = Vec3(0.0);
    pic = p.velocity;

    WeightOverParticleNeighbourhood(
        mParams, p,
        [&](IVec3 pos, Float weight) {
            // Transfer mass
            const Cell& cell = g.Get(pos.x, pos.y, pos.z);
            flip += cell.VelocityStar * weight;
            pic += (cell.VelocityStar - cell.Velocity) * weight;
        }
    );
}

void ParticleSystem::UpdateVelocities(const Grid& g, MTIterator& mt)
{
    mt.IterateOverVector(mParticles, [&](Particle& p) {
        Vec3 flip;
        Vec3 pic;

        CalculateFlipPicVelocity(p, g, pic, flip);
        p.velocity = (1 - mParams.ALPHA) * pic + mParams.ALPHA * flip;
    });
}


void ParticleSystem::UpdatePositions(Float dt, MTIterator& mt) 
{
    mt.IterateOverVector(mParticles, [&](Particle& p) {
        p.pos += dt * p.velocity; 
    });
}

std::vector<Particle>& ParticleSystem::GetParticles()
{
    return mParticles;
}

const std::vector<Particle>& ParticleSystem::GetParticles() const
{
    return mParticles;
}

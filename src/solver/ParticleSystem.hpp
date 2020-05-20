#pragma once

#include <vector>

#include "Common.hpp"
#include <array>
#include "SimulationParameters.hpp"

using ParticleHandle = Uint;

class Grid;
class ParticleSystem;
class MTIterator;

class Particle {
public:
    Vec3 pos;
    Float mass;
    Vec3 velocity;
    Float volume;
    Mat3 m_F_p;
    Mat3 m_F_e;
    Mat3 m_R_e;
    
    // We cache this since it ends up being quite expensive
    // to do this every time (calculating the gradient is 27 branches )
    Uint numNeighbours;
    std::array<IVec3, 9> neighbours_coords;
    std::array<Float, 9> neighbours_nx;
    std::array<Vec3, 9> neighbours_nxgrad;

private:
    Particle(const Vec3& pos, Float mass, const Vec3& velocity);
    friend ParticleSystem;
};

class ParticleSystem {
public:
    ParticleSystem(const SimulationParameters& parameters);
    void AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass);

    Mat3 CalculateCauchyStress(const Particle& p) const;

    void CacheParticleGrads(const Grid& g, MTIterator& mt);
    void EstimateParticleVolumes(const Grid& g, MTIterator& mt);
    void UpdateDeformationGradients(const Float dt, const Grid& g, MTIterator& mt);
    void UpdateVelocities(const Grid& g, MTIterator& mt);
    void BodyCollisions(Float dt, MTIterator& mt);
    void UpdatePositions(Float dt, MTIterator& mt);

    // todo:  This interface is 'dirty' - does a better way for contignuous particle data access exist?
    std::vector<Particle>& GetParticles();
    const std::vector<Particle>& GetParticles() const;

private:
    void CalculateFlipPicVelocity(const Particle& p, const Grid& g, Vec3& flip, Vec3& pic) const;
    Mat3 CalculateVelocityGradient(const Particle& p, const Grid& g) const;

    const SimulationParameters& mParams;
    std::vector<Particle> mParticles;
};

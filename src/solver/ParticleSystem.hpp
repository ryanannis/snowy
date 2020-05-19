#pragma once

#include <vector>

#include "Common.hpp"
#include "SimulationParameters.hpp"

using ParticleHandle = Uint;

class Grid;
class ParticleSystem;

class Particle {
public:
    Vec3 pos;
    Float mass;
    Vec3 velocity;
    Float volume;
    Mat3 m_F_p;
    Mat3 m_F_e;
    Mat3 m_R_e;

private:
    Particle(const Vec3& pos, Float mass, const Vec3& velocity);
    friend ParticleSystem;
};

class ParticleSystem {
public:
    ParticleSystem(const SimulationParameters& parameters);
    void AddParticle(const Vec3& pos, const Vec3& velocity, const Float mass);

    Mat3 CalculateCauchyStress(const Particle& p) const;
    
    void EstimateParticleVolumes(Grid& g);
    void UpdateDeformationGradients(const Float dt, const Grid& g);
    void UpdateVelocities(const Grid& g);
    void BodyCollisions(Float dt);
    void UpdatePositions(Float dt);

    // todo:  This interface is 'dirty' - does a better way for contignuous particle data access exist?
    std::vector<Particle>& GetParticles();
    const std::vector<Particle>& GetParticles() const;

private:
    Vec3 CalculatePICVelocity(const Particle& p, const Grid& g) const;
    void CalculateFlipPicVelocity(const Particle& p, const Grid& g, Vec3& flip, Vec3& pic, bool fag) const;
    Mat3 CalculateVelocityGradient(const Particle& p, const Grid& g) const;

    const SimulationParameters& mParams;
    std::vector<Particle> mParticles;
};

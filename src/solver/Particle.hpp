#pragma once

#include "Common.hpp"

class Grid;

class Particle {
public:
    Particle(const Vec3& pos);
    Mat3 CalculateCauchyStress() const;
    void UpdateDeformationGradient(Float dt, Vec3 dv, const Grid& g);
    void UpdateVelocity(const Grid& g);
    void UpdatePosition(Float dt);

    Vec3 Position;
    Float Mass = 0;
    Vec3 Velocity = {};
    Float Volume = 0;

private:
    Vec3 CalculatePICVelocity(const Grid& g) const;
    void CalculateFlipPicVelocity(const Grid& g, Vec3& flip, Vec3& pic) const;
    Mat3 CalculateVelocityGradient(const Grid& g) const;

    Mat3 m_F_e;
    Mat3 m_R_e;
    Mat3 m_F_p;
};

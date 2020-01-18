#include "Particle.hpp"

#include "Grid.hpp"
#include "Math.hpp"
#include "glm/gtx/matrix_operation.hpp"

Particle::Particle(const Vec3& pos) :
    Position(pos)
{
}

Mat3 Particle::CalculateVelocityGradient(const Grid& g) const 
{
    Mat3 velGrad = {};
    assert(velGrad == Mat3(Float(0.0))); // remove me

    for (int i = -1; i <= 2; i++)
    {
        for (int j = -1; j <= 2; j++)
        {
            for (int k = -1; k <= 2; k++)
            {
                // The coordinate of the cell we are rasterizing to
                int ix = static_cast<int>(Position.x) + i;
                int iy = static_cast<int>(Position.y) + j;
                int iz = static_cast<int>(Position.z) + k;

                Vec3 nxGrad = gridWeightGrad(Position.x, ix, Position.y, iy, Position.z, iz);

                velGrad += glm::outerProduct(g.Get(ix, iy, iz).VelocityStar, nxGrad); 
            }
        }
    }

    return velGrad;
}

Mat3 Particle::CalculateCauchyStress() const
{
    Float j_p = glm::determinant(m_F_p);
    Float j_e = glm::determinant(m_F_e);

    Float hardening = exp(HARDENING * (1 - j_p));
    Float mu = MU_0 * hardening;
    Float lambda = LAMBDA_0 * hardening;
    return Float(2.0) * mu * (m_F_e - m_R_e) * glm::transpose(m_F_e) + Mat3(lambda * (j_e - 1) * j_e);
}
 

void Particle::UpdateDeformationGradient(Float dt, Vec3 dv, const Grid& g)
{
    // First attribute all new changes to elastic part of deformation
    auto new_fe = (Mat3(Float(1.0)) + dt * CalculateVelocityGradient(g)) * m_F_e;
    
    Mat3 u;
    Mat3 s;
    Mat3 sinv(Float(0.0));
    Mat3 v;
    svd3(new_fe, u, s, v);

    for(Uint i = 0; i < 3; i++) {
        s[i][i] = glm::clamp(s[i][i], Float(1.0) - PHI_C, Float(1.0) + PHI_S);
        sinv[i][i] = Float(1.0) / s[i][i];
    }

    m_F_e = u * s * glm::transpose(v);
    m_F_p = v * sinv * glm::transpose(u) * m_F_p;
}

void Particle::CalculateFlipPicVelocity(const Grid& g, Vec3& flip, Vec3& pic) const
{
    flip = Vec3(0.0);
    pic = Vec3(0.0);

    for (int i = -1; i <= 2; i++)
    {
        for (int j = -1; j <= 2; j++)
        {
            for (int k = -1; k <= 2; k++)
            {
                // The coordinate of the cell we are rasterizing to
                int ix = static_cast<int>(Position.x) + i;
                int iy = static_cast<int>(Position.y) + j;
                int iz = static_cast<int>(Position.z) + k;

                const Float nx = gridWeight(Position.x, ix, Position.y, iy, Position.z, iz);
                const Cell& cell = g.Get(ix, iy, iz);

                flip += cell.VelocityStar * nx;
                pic += (cell.VelocityStar - cell.Velocity) * nx;
            }
        }
    }
}

void Particle::UpdateVelocity(const Grid& g)
{
    Vec3 flip;
    Vec3 pic;

    CalculateFlipPicVelocity(g, flip, pic);
    Velocity += (1 - ALPHA) * pic + ALPHA * flip;
}


void Particle::UpdatePosition(Float dt)
{
    Position += dt * Velocity; 
}
#pragma once

#include <cstdint>
#include <glm/glm.hpp>

// Numeric types
using Uint = uint64_t;
using Float = float;
using Vec3 = glm::vec3;
using Mat3 = glm::mat3;

static const Float H = 1.0; // cell size
static const Float HARDENING = 1.0;
static const Float MU_0 = 1.0;
static const Float LAMBDA_0 = 1.0;
static const Float PHI_C = 0.025;
static const Float PHI_S = 0.0075;
static const Float ALPHA = 0.95;

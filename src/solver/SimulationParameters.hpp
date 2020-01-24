#pragma once

#include "Common.hpp"

struct SimulationParameters {
    Float H = 0.001; // cell size
    Float HARDENING = 1.0;
    Float MU_0 = 1.0;
    Float LAMBDA_0 = 1.0;
    Float PHI_C = 0.025;
    Float PHI_S = 0.0075;
    Float ALPHA = 0.95;
    Float GRAVITY = -9.81;
};

#pragma once

#include "Common.hpp"

struct SimulationParameters {
    Float H = 1.0; // cell size
    Float HARDENING = 10.0;
    Float MU_0 = 58333.0;
    Float LAMBDA_0 = 38888.0;
    Float PHI_C = 0.025;
    Float PHI_S = 0.005;
    Float ALPHA = 0.95;
    Float GRAVITY = -9.81;
};

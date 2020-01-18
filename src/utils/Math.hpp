#pragma once

#include <cstdint>
#include "Common.hpp"

void svd3 (const Mat3 mat, Mat3& u, Mat3& s, Mat3& v);

// Numeric types
static const double EPSILON = 0.000000000001;


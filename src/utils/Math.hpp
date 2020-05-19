#pragma once

#include <cstdint>
#include "Common.hpp"

// todo:  "simlation math" should be seperated out from "general math"

void svd3(const Mat3& mat, Mat3& u, Mat3& s, Mat3& v);
Vec3 gridWeightGrad(Float h, Float x, Float ix, Float y, Float iy, Float z, Float iz);
Float gridWeight(Float h, Float x, Float ix, Float y, Float iy, Float z, Float iz);

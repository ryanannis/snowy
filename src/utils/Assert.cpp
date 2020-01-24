#include "assert.hpp"

void ASSERT_VALID_FLOAT(Float f)
{
	assert(std::isfinite(f));
}

void ASSERT_VALID_VEC3(Vec3 v)
{
	assert(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
}

void ASSERT_VALID_MAT3(Mat3 v)
{
	assert(
		std::isfinite(v[0][0]) && std::isfinite(v[0][1]) && std::isfinite(v[0][2]) &&
		std::isfinite(v[1][0]) && std::isfinite(v[1][1]) && std::isfinite(v[1][2]) &&
		std::isfinite(v[2][0]) && std::isfinite(v[2][1]) && std::isfinite(v[2][2])
	);
}

void ASSERT_NOT_SINGULAR(Mat3 v)
{
}

#include "gtest/gtest.h"
#include "CPUSolver.hpp"

#include <iostream>

TEST(rasterization, roundtrip)
{
    RecordProperty("description", "Makes sure a roundtrip of rasterizing and derasterizing produces the same result.");
    Grid g(10, 10, 10);
    Particle p(glm::vec3(3.5, 3.5, 3.5));
    p.Mass() = 1.0;
    RasterizationUtils::RasterizeMassToGrid(g, { p });
    std::cout << g << std::endl;
};
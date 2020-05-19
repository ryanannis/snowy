#include "gtest/gtest.h"
#include "CPUSolver.hpp"
#include "Math.hpp"

#include <iostream>

// Make sure amount of energy remains constant when transferred to grid and back
TEST(IntegrationTests, ConstantEnergy) {
    Mat3 u(1);
    Mat3 v(1);
    Mat3 s(50);
    Mat3 base(
        2.0, 1.0, 0.0,
        1.0, 2.0, 0.0,
        2.0, 0.0, 3.0
    );

    svd3(base,u,s,v);

    Mat3 f = u * s * glm::transpose(v);
    Mat3 matdiff = base - f;


    std::cout << matdiff[0][0] << " " << matdiff[1][0] << " " << matdiff[2][0] << std::endl
              << matdiff[0][1] << " " << matdiff[1][1] << " " << matdiff[2][1] << std::endl
              << matdiff[0][2] << " " << matdiff[1][2] << " " << matdiff[2][2] << std::endl;
     
}
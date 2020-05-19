#include "Common.hpp"

#include <string>
#include <vector>

struct Triangle
{ 
    Vec3 v1;
    Vec3 v2;
    Vec3 v3;
};

class Mesh
{
public:
    Mesh(const std::string& path);

private:
    void initFromObj();

    std::vector<Triangle> mTriangle;
};


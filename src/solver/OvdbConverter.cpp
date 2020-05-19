#define NOMINMAX 1

#include "OvdbConverter.hpp"
#include "Math.hpp"

#include <iostream>
#include <vector>
#include <openvdb/openvdb.h>

// Private helper functions

// Implementation

OvdbConverter::OvdbConverter(std::shared_ptr<SimulationOutput> data) :
    mData(data)
{
}

void OvdbConverter::Output(const std::string& path) const
{
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);
    const Float H = 1;

    grid->setTransform(openvdb::math::Transform::createLinearTransform(H)); // TODO:  This should be H

    // Identify the grid as a level set.
    // todo:  why is the examples using the level set for smoke... ?  can that somehow be represented as a narrow band level set?
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    // Name the grid "LevelSetSphere".
    grid->setName("snowdensity");

    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    for(const Particle &p : mData->GetParticles())
    {
        const Vec3& Position = p.pos;
        // Use our rasterization kernel to accumulate density onto the grid

        // todo:  this is a nasty thing to copy and paste
        // i'm guessing that the grid random access is awful and it's probably better to use the normal accumulator and
        // then iterating that (considering I don't have sparsity implemented in the first place)
        
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                for (int k = -2; k < 3; k++)
                {
                    // The coordinate of the cell we are rasterizing to
                    int ix = static_cast<int>(Position.x / H) + i;
                    int iy = static_cast<int>(Position.y / H) + j;
                    int iz = static_cast<int>(Position.z / H) + k;

                    Float nx = gridWeight(H, Position.x / H, ix, Position.y / H, iy, Position.z / H, iz);

                    openvdb::Coord xyz(ix, iy, iz);
                    accessor.setValue(xyz, accessor.getValue(xyz) + nx * p.mass);
                }
            }
        }
    }


    // Create a VDB file object.
    openvdb::io::File file(path);

    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
}
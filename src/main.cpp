#include <iostream>
#include <fstream>
#include <string>

#include "CPUSolver.hpp"

// Our interfacing is done through file (currently)
int main(int argc, char* argv[]) {
    SimulationParameters sp;
    sp.H = 1.0;
    
    CPUSolver solver(
        IVec3(50, 50, 50),
        1.0 / 24.0,
        sp
    );

    // Add a cube of snow
    Vec3 center = Vec3(25, 25, 25);
    for (int x = -5; x < 5; x++) {
        for (int y = -5; y < 5; y++) {
            for (int z = -5; z < 5; z++) {
                solver.AddParticle(
                    center + Vec3(Float(x), Float(y), Float(z)),
                    Vec3(Float(0.0), Float(0.0), Float(0.0)),
                    0.001 // 1g per particle
                );
            }
        }
    }

    for (int i = 0; i < 24; i++)
    {
        std::cout << "outputting frame " + std::to_string(i) << std::endl;
        std::ofstream outfile("out" + std::to_string(i) + ".txt");
        for (Particle p : solver.NextFrame()) {
            outfile << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << std::endl;
        }
        outfile.close();
    }

    return 0;
}
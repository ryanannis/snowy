#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>

#include "CPUSolver.hpp"
#include "OvdbConverter.hpp"

// Our interfacing is done through file (currently)
int main(int argc, char* argv[]) {
    std::srand(std::time(nullptr));
    SimulationParameters sp;
    sp.H = 1.0; // 1m grid size 
    
    CPUSolver solver(
        IVec3(125, 60, 60),
        1.0 / 24.0,
        sp
    );

    // Add a cube of snow
    Vec3 center = Vec3(10, 30, 30);
    for (int x = -30; x < 30; x++) {
        for (int y = -30; y < 30; y++) {
            for (int z = -30; z < 30; z++) {
                if (Float(x) * Float(x) + Float(y) * Float(y) + Float(z) * Float(z) > 9.0) {
                    continue;
                }

                float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                float r3 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                
                solver.AddParticle(
                    center + Vec3(Float(x) / 6, Float(y) / 6.0, Float(z) /6.0),
                    Vec3(Float(250.0), Float(0.0), Float(0.0)),
                    200.0 // 1g per particle
                );
            }
        }
    }

    center = Vec3(35, 30, 30);
    for (int x = -3; x < 0; x++) {
        for (int y = -90; y < 90; y++) {
            for (int z = -90; z < 90; z++) {
                
                solver.AddParticle(
                    center + Vec3(Float(x) / 6.0, Float(y) / 6.0, Float(z) / 6.0),
                    Vec3(Float(-15.0), Float(1.0), Float(-0.5)),
                    1.0 // 1g per particle
                );
            }
        }
    }

    std::shared_ptr<SimulationOutput> simoutput;

    for (int i = 0; i < 48; i++)
    {
        std::cout << "outputting frame " + std::to_string(i) << std::endl;
        std::ofstream outfile("out" + std::to_string(i) + ".txt");
        solver.NextFrame();
        simoutput = solver.GetOutput();
        for (Particle p : simoutput->GetParticles()) {
            outfile << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << std::endl;
        }
        OvdbConverter converter(simoutput);
        converter.Output("sim_" + std::to_string(i) + ".vdb");

        outfile.close();
    }

    return 0;
}
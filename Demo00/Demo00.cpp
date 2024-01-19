
// Original nBody Demo without optimizations

#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

struct Particle {
    double x, y, z;
    double mass;
    double vx, vy, vz;
};

class ParticleSystem {
public:
    std::vector<Particle> particles;

    ParticleSystem(int num_particles) {
        for (int i = 0; i < num_particles; i++) {
            Particle p = { rand() % 1000, rand() % 1000, rand() % 1000,
                          1.0,
                          0.0, 0.0, 0.0 };
            particles.push_back(p);
        }
    }

    void update() {
        for (int i = 0; i < particles.size(); i++) {
            for (int j = i + 1; j < particles.size(); j++) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;

                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                double force = particles[i].mass * particles[j].mass / (dist * dist * dist);

                particles[i].vx += force * dx;
                particles[i].vy += force * dy;
                particles[i].vz += force * dz;

            	particles[j].vx -= force * dx;
                particles[j].vy -= force * dy;
                particles[j].vz -= force * dz;
            }
        }

        for (auto& p : particles) {
            p.x += p.vx;
            p.y += p.vy;
            p.z += p.vz;
        }
    }
};



int main(int argc, char* argv[])
{
    int numberOfParticles = 1000;
    int numberOfSteps = 100;

    if (argc > 1) {
        numberOfParticles = std::atoi(argv[1]);
        if (argc == 3) {
            numberOfSteps = std::atoi(argv[2]);
        }
    }

    ParticleSystem system(numberOfParticles);

    // start timing
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numberOfSteps; i++) {
        system.update();
    }

    // stop timing
    auto stop = std::chrono::high_resolution_clock::now();

    // compute the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    //std::cout << std::endl << "Runtime: " << duration.count() << " microseconds" << std::endl;
    std::cout << "\n\tOriginal nBody Demo without optimizations...\n";
    std::cout << "\t\tNumParticles\t\tNumberOfSteps\t\tRuntime(microseconds)\n";
    std::cout << "\t\t" << numberOfParticles << "\t\t\t" << numberOfSteps << "\t\t\t" << duration.count();

    return 0;
}

#include <iostream>
#include <vector>
#include <cmath>
#include <execution>
#include <random>
#include <array>

class ParticleSystem
{
public:
    std::vector<std::array<double, 3>> positions;
    std::vector<std::array<double, 3>> velocities;
    std::vector<std::array<double, 3>> accelerations;
    std::vector<double> masses;

    ParticleSystem(int num_particles)
	{
        // Reserve space for all particles up front
        positions.reserve(num_particles);
        velocities.reserve(num_particles);
        accelerations.reserve(num_particles);
        masses.reserve(num_particles);

        // Use a better random number generation scheme
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1000.0);

        // Initialize particles...
        for (int i = 0; i < num_particles; ++i) {
            std::array<double, 3> pos = { {1e18 * exp(-1.8) * (0.5 - dis(gen)),
                                          1e18 * exp(-1.8) * (0.5 - dis(gen)),
                                          1e18 * exp(-1.8) * (0.5 - dis(gen))} };

            std::array<double, 3> vel = { {1e18 * exp(-1.8) * (0.5 - dis(gen)),
                                          1e18 * exp(-1.8) * (0.5 - dis(gen)),
                                          1e18 * exp(-1.8) * (0.5 - dis(gen))} };

            double mass = 1.98892e30 * dis(gen) * 10 + 1e20;

            positions.push_back(pos);
            velocities.push_back(vel);
            accelerations.push_back({ 0.0, 0.0, 0.0 });
            masses.push_back(mass);
        }
    }

    void update()
	{
        // Reset accelerations to zero
        for (auto& acc : accelerations) {
            acc.fill(0.0);
        }

        // Calculate forces in parallel using C++20 parallel algorithms
        std::for_each(std::execution::par, positions.begin(), positions.end(),
            [&](std::array<double, 3>& pi_pos, size_t i) {
                for (size_t j = 0; j < positions.size(); ++j) {
                    if (i != j) {
                        std::array<double, 3>& pj_pos = positions[j];
                        double dx = pj_pos[0] - pi_pos[0];
                        double dy = pj_pos[1] - pi_pos[1];
                        double dz = pj_pos[2] - pi_pos[2];

                        double dist_squared = dx * dx + dy * dy + dz * dz;
                        double dist_sixth = dist_squared * dist_squared * dist_squared;
                        double force = masses[i] * masses[j] / (dist_sixth + 1e-12);

                        accelerations[i][0] += force * dx;
                        accelerations[i][1] += force * dy;
                        accelerations[i][2] += force * dz;
                    }
                }
            }
        );

        // Update velocity and position
        for (size_t i = 0; i < positions.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                velocities[i][j] += accelerations[i][j];
                positions[i][j] += velocities[i][j];
            }
        }
    }
};


int main(int argc, char** argv)
{
    int numberOfParticles = 1000;
    int numberOfSteps = 100;

    // To do... better validation of args.
    if (argc > 1) {
        numberOfParticles = std::atoi(argv[1]);
        if (argc == 3) {
            numberOfSteps = std::atoi(argv[2]);
        }
    }

    try {
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

        std::cout << "\n\tModified nBody Demo with C++20 using an 'Array of Structures' to a 'Structure of Arrays' for Particle data...\n";
        std::cout << "\t\tNumParticles\t\tNumberOfSteps\t\tRuntime(microseconds)\n";
        std::cout << "\t\t" << numberOfParticles << "\t\t\t" << numberOfSteps << "\t\t\t" << duration.count();
    }
    catch (std::exception& e)
    {
        std::cout << "Exception: " << e.what() << "\n";
        std::terminate();
    }
    return 0;
}




#include <iostream>
#include <vector>
#include <cmath>
#include <execution>
#include <random>
#include <cstdlib>


// 3D vector class

typedef struct Particle
{
    Particle() { init(); }

    void init() {
        pos[0] = 0.; pos[1] = 0.; pos[2] = 0.;
        vel[0] = 0.; vel[1] = 0.; vel[2] = 0.;
        acc[0] = 0.; acc[1] = 0.; acc[2] = 0.;
        mass = 0.;
    }

    double pos[3];        // position
    double vel[3];		  // velocity
    double acc[3];		  // acceleration
    double mass;
} Particle;


class ParticleSystem
{
public:
    std::vector<Particle> particles;

    ParticleSystem(int num_particles)
	{
        //
		// Optimization:  reserve space for all particles up front.
		//
        particles.reserve(num_particles);

        // Use a better random number generation scheme than rand()	
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1000.0);

        // Initialize particles...
    	for (int i = 0; i < num_particles; ++i) {
    		Particle p;
        	// random position
			p.pos[0] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			p.pos[1] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			p.pos[2] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			// random velocity
			p.vel[0] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			p.vel[1] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			p.vel[2] = 1e18 * exp(-1.8) * (0.5 - dis(gen));
			// random mass
			p.mass = 1.98892e30 * dis(gen) * 10 + 1e20;
			// initial force
			p.acc[0] = 0.0;
			p.acc[1] = 0.0;
			p.acc[2] = 0.0;
            particles.push_back(p);
		}
    }

    void update()
	{
        // Reset accelerations to zero
        for (auto& p : particles) {
            p.acc[0] = p.acc[1] = p.acc[2] = 0.0;
        }

        // Calculate forces in parallel using C++20 parallel algorithms
        std::for_each(std::execution::par, particles.begin(), particles.end(),
            [&](Particle& pi) {
                for (auto& pj : particles) {
                    if (&pi != &pj) {
                        double dx = pj.pos[0] - pi.pos[0];
                        double dy = pj.pos[1] - pi.pos[1];
                        double dz = pj.pos[2] - pi.pos[2];

                        //
                        // Optimization: Avoid redundant calls to std::sqrt in this inner loop
                        //
                        double dist_squared = dx * dx + dy * dy + dz * dz;
                        double dist_sixth = dist_squared * dist_squared * dist_squared;
                        double force = pi.mass * pj.mass / (dist_sixth + 1e-12); // Small epsilon to avoid division by zero

                        pi.acc[0] += force * dx;
                        pi.acc[1] += force * dy;
                        pi.acc[2] += force * dz;
                    }
                }
            }
        );

        // Update velocity and position
        for (auto& p : particles) {
            for (int i = 0; i < 3; ++i) {
                p.vel[i] += p.acc[i];
                p.pos[i] += p.vel[i];
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

        std::cout << "\n\tModified nBody Demo with C++20 standard library optimizations...\n";
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

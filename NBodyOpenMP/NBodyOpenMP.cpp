

/*
 *  NBodyOpenMP.cpp
 *
 *  Demonstration of OpenMP library parallel algorithms applied to the n-body problem.
 *
 */

#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>

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

		//
		//  Thoughts on thread safety...
		//
		//  The following code is not using a Parallel for loop because the std::vector::push_back() is not thread-safe.
		//  And, the std::random_device is not thread-safe.  We could rewrite this so that each thread has its own
		//  std::random_device, but that would be inefficient too.  Research further...
		//

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
		int n = particles.size();

		// Reset accelerations to zero
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			particles[i].acc[0] = particles[i].acc[1] = particles[i].acc[2] = 0.0;
		}

		// Calculate forces in parallel using OpenMP
#pragma omp parallel for 
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (i != j) {
					double dx = particles[j].pos[0] - particles[i].pos[0];
					double dy = particles[j].pos[1] - particles[i].pos[1];
					double dz = particles[j].pos[2] - particles[i].pos[2];

					double dist_squared = dx * dx + dy * dy + dz * dz;
					double dist_sixth = dist_squared * dist_squared * dist_squared;
					double force = particles[i].mass * particles[j].mass / (dist_sixth + 1e-12);

					particles[i].acc[0] += force * dx;
					particles[i].acc[1] += force * dy;
					particles[i].acc[2] += force * dz;
				}
			}
		}

		// Update velocity and position
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < 3; ++k) {
				particles[i].vel[k] += particles[i].acc[k];
				particles[i].pos[k] += particles[i].vel[k];
			}
		}
	}
};


int main(int argc, char** argv)
{
	int numberOfParticles = 10000;
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

		std::cout << "\n\tModified nBody Demo with OpenMP library optimizations...\n";
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


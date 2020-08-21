#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <assert.h>

using namespace std;

SPH_main domain;

// Test that velocity magnitude is lower than the speed of sound
void test_v(SPH_particle* part, SPH_main* domain)
{
	double v_mag = sqrt(part->v[0] * part->v[0] + part->v[1] * part->v[1]);
	assert(v_mag < domain->c0);
}

// Test that density value is reasonable
void test_rho(SPH_particle* part, SPH_main* domain)
{
	assert(part->rho < domain->rho0 + 500 && part->rho > domain->rho0 - 500);
}

// Test that particle is not far outside simulation domain
void test_position(SPH_particle* part, SPH_main* domain)
{
	assert(part->x[0] > domain->min_x[0] - (4 * domain->h) && part->x[0] < domain->max_x[0] + (2 * domain->h));
	assert(part->x[1] > domain->min_x[1] - (4 * domain->h) && part->x[1] < domain->max_x[1] + (2 * domain->h));
}

int main(void)
{
	string geometry = "default";

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid
	domain.place_points(domain.min_x, domain.max_x, geometry);			//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

	// Initialising timestepping
	double dt = 0.1 * domain.h / domain.c0; // first timestep
	double tmax = 3.0; // total simulation time
	double time = 0.; // current time

	int update_freq = 10; // smoothing density
	//double file_time = 0.; // times to output files, gets updated in loop

	// Start timer
	auto start = std::chrono::system_clock::now();

	int step = 0;
	//int count = 0;
	while (time <= tmax)
	{
		domain.allocate_to_grid();

		// Smooth density
		if (step % update_freq == 0)
		{
			#pragma omp parallel for schedule (guided)		// Parallelize the loop
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.smooth(&domain.particle_list[i], dt);
		}

		// Iterate through neighbors and fill in information for time-stepping, for every particle
		#pragma omp parallel for schedule (guided)		// Parallelize the loop
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.neighbour_iterate(&domain.particle_list[i], dt);

		// Update values of every particle
		#pragma omp parallel for schedule (guided)		// Parallelize the loop
		for (int i = 0; i < domain.particle_list.size(); i++)
			//domain.forward_euler(&domain.particle_list[i], dt);
			domain.predictor_corrector_half(&domain.particle_list[i], dt);

		domain.allocate_to_grid();

		// Iterate through neighbors and fill in information for time-stepping, for every particle
		#pragma omp parallel for schedule (guided)		// Parallelize the loop
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.neighbour_iterate(&domain.particle_list[i], dt);

		// Update values of every particle
		#pragma omp parallel for schedule (guided)		// Parallelize the loop
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.predictor_corrector_full(&domain.particle_list[i], dt);
			test_position(&domain.particle_list[i], &domain);
			test_v(&domain.particle_list[i], &domain);
			test_rho(&domain.particle_list[i], &domain);
		}

		time += dt;
		step++;

		// Output VTP file every 0.1s
		/*if (time >= file_time)
		{
			string fname = "./tests/bin/example_" + to_string(count) + ".vtp";
			write_file(fname.c_str(), &domain.particle_list);
			count++;
			file_time += 0.1;
		}*/

		// Update new timestep
		domain.get_dt();
		dt = domain.dt;
	}
	auto end = std::chrono::system_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end - start);

	if (duration.count() > tmax * 100)
	{
		cerr << "\n--------------------------------------------------------------------\n\n";
		cerr << "Runtime has exceeded " << tmax * 100 << " seconds. Simulation is unstable!";
		cerr << "\n\n--------------------------------------------------------------------\n";
	}

	cout << "\n\nTime taken for simulating " << tmax << " s: " << duration.count() << " seconds" << endl;
	cout << "\nAll tests passed.\n";
	return 0;
}

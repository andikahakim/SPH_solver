#pragma once

#include <vector>
#include "SPH_2D.h"

/// Function to write simulation result to .vtp
int write_file(const char* filename,
	std::vector<SPH_particle>* particle_list);

# ACSE-4 Project 3: Smooth Particle Hydrodynamics Solver
#### Team Katrine: Andika Hakim, Qingyang Li, Helen Situ, Chao-Lun Liu, Chen Qian, Nandong Yan

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

## Animation

[![Youtube simulation video]](https://www.youtube.com/watch?v=POnmzzhc5E0)

[![Watch the video]](https://github.com/acse-2019/acse-4-sph-katrine/blob/master/video/animation.mp4)

## Compilation/Installation Guide

To install the module `acse-4-sph-katrine
` clone the respository to your computer with 
``
git clone https://github.com/acse-2019/acse-4-sph-katrine.git
``
Then change directory into the root of the repo you just cloned 
``
cd acse-4-sph-katrine
``
and the run 
```
make
```
The make command will automatically fetch the ``Makefile`` instructions in the data and execute the following scripts
```
g++ -o build/SPH_2D.o -c src/SPH_2D.cpp -Wall -std=c++17 -fopenmp -fpermissive -Iincludes
g++ -o build/SPH_Snippet.o -c src/SPH_PC_Snippet.cpp -Wall -std=c++17 -fopenmp -fpermissive -Iincludes
g++ -o build/file_writer.o -c src/file_writer.cpp -Wall -std=c++17 -fopenmp -fpermissive -Iincludes
g++ -o bin/SPH_2D build/SPH_2D.o build/SPH_Snippet.o build/file_writer.o
```
**NOTE: By default, the Predictor-Corrector scheme snipped will be compiled. However, the alternate geometries are only implemented and fully functional in the Forward Euler scheme. To compile the Forward Euler scheme snippet, please rename the alternate `FE_makefile` to `makefile` and replace the current `makefile`.

If users want to run the repo on command line, you can run below instruction, example:

```
g++ -Xpreprocessor -fopenmp -lomp -o main SPH_FE_Snippet.cpp SPH_2D.cpp file_writer.cpp -std=c++17

./main
```

## User instructions

Save
in the file ``/src`` and  ``/icnludes``

-- [SPH_Snippet.cpp](SPH_2D.h)
-- [SPH_FE_Snippet.cpp](SPH_2D.h)
-- [SPH_PC_Snippet.cpp](SPH_2D.h)
-- [SPH_2D.cpp](SPH_2D.h)
-- [file_writer.cpp](file_writer.h)

to your code directory and inform it as a header in your file.

```C++
#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iomanip>
#include <cmath>
#include <complex>
#include <math.h>
#include <fstream>
#include <vector>
#include <omp.h>
#define _USE_MATH_DEFINES
```

When dealing with a main, which will call write_file function:

```C++
#include "file_writer.h"
#include <string>
#include <fstream>
#include <vector>
```

Developer users with experience in the fields of C ++ and SPH will be able to take our code and build on top of it. For these people, we encourage you to read through the code and understand how it works, but to get you started, you can use the following lines of code in a C ++ script to initialize and run your first simulation:

### OpenMP  : 

OpenMP consists of a set of compiler #pragmas that control how the program works. The pragmas are designed so that even if the compiler does not support them, the program will still yield correct behavior, but without any parallelism.

C++

STL is not thread-safe. If you use STL containers in a parallel context, you must exclude concurrent access using locks or other mechanisms. Const-access is usually fine, as long as non-const access does not occur at the same time.

Exceptions may not be thrown and caught across omp constructs. That is, if a code inside an ``omp for`` throws an exception, the exception must be caught before the end of the loop iteration; and an exception thrown inside a parallel section must be caught by the same thread before the end of the ``parallel``section.

GCC

``fork()`` is troublematic when used together with OpenMP. See the chapter "OpenMP and fork()" above for details.

### Optional arguments  :

```C++
/* Initialising timestepping */
double dt = 0.1 * domain.h / domain.c0; // first timestep
double tmax = 30.0; // total simulation time
double time = 0.; // current time

int update_freq = 10; // smoothing density
double file_time = 0.; // times to output files, gets updated in loop
```

### Create boundary & Set the grid points :

```C++

SPH_main domain;

//Set simulation parameters
domain.set_values();
//initialise simulation grid
domain.initialise_grid();
//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
domain.place_points(domain.min_x, domain.max_x, "Bubble");
```

## Documentation

```

Information on the implementation and functions included are in the /documentation/html ``index.html`` document.

```
## Structure

### The SPH_2D module:

The sph_2D module is mainly used for generating a  domain and dynamic particles (and their properties) in the meshed domain.

 ### The  SPH_main class: 
 The SPH_main class, within "SPH_2D.cpp" and "SPH_2D.h" used to initialise a single object of the main SPH type (a initial domain containing boundaries and liquid particles) <br>

#### Class methods:

##### **SPH_main():**

- This is the class object that creates the domain onto which the particles will be placed

**set_values():**

- This method sets the simulation domain parameters.
- Domain reference values such as the fluid speed of sound, the domain reference density and the solution end time are given .<br>

**initialise_grid():**

- Initialises a simulation grid and calculates the size of the array needed for an attribute called search_grid[]. search_grid[] is used later on to calculate the contributions of neighbouring particles in neighbour_iterate().

**place_points(min,  max, geometry):**

- within the initialised grid, a particle at every increment of dx is created using the SPH_particle() class. The particles are assigned an index relative to the search grid and then appended into a list called particle_list[].

**allocate_to_grid()**

- particles are assigned to a search grid using the index created earlier.

**neighbour_iterate(part, dt)**

* this function iterates through every particle in the domain. It quantifies the contributions of the neighbours by using it's search_grid index and then approximates the derivatives in the navier-stokes equation using a smoothing kernel

**forward_euler(part, dt)**

-  this function time steps the particles in the particle list. At every time increment, neighbour_iterate and update_values are used to update the particle parameters for the next time increment
-  at the end of the forward_wrapper(), the each particle_list at prespecified time increment is combined and then saved into a file called 'State.npy'.

**predictor_corrector_half(part, dt)**
**predictor_corrector_full(part, dt)**

- The above functions allow a particle’s acceleration and rate of change of density to be calculated as a function of its neighbour’s positions, velocities and densities/pressures. We now need to be able to update the particle’s position, velocity and density, with the updated density being used to calculate a new pressure. To save complexity we will be doing this explicitly.

- The added complexity with the predictor-corrector scheme is that 2 values for each of the variables need to be stored (the initial value and the half-step value) and the neighbour searching needs to be based on the appropriate 𝐱 values.


### The SPH_particle Class:

The SPH_particle class, within "SPH_2D.h" is a class which creates an object of a single particle with attributes. The class contains the attributes including position (x), velocity (v), acceleration (a), density (rho), pressure (P) and the instantaneous change in density (D).

Within the SPH_particle class are class functions:

**calc_index()**

- this function gives assigns a list number to the current particle. This list number is used to allocate the particle to the search grid


## Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
python run_tests.py
```

Testing is performed using ``run_tests.py``.  Four sanity checks are implemented:

-  test_v() : Test that velocity magnitude is lower than the speed of sound

- test_rho() : Test that density value is reasonable

- test_position() : Test that particle is not far outside simulation domain


## Contributors

- [**Professor Stephen Neethling**](https://www.imperial.ac.uk/people/s.neethling)
- [**Kartrine group**](https://github.com/acse-2019/acse-4-sph-katrine)

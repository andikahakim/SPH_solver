#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class SPH_main;

/// Class object that contain particle information
///
/// Every particle on the system is having unique information that is stored in this class
class SPH_particle
{
public:
    /// Constructor for particle class
    SPH_particle();

    /// Calculating grid allocation index
    ///
    /// Updates list_num, an array that stores the grid row and grid column indices to allocate particle to grid cells.
    /// @return void
    void calc_index(void);


    /// Variable to store position (vector - x and y direction) at one time step
    double x[2];
    /// Variable to store velocity (vector - x and y direction) at one time step
    double v[2];
    /// Variable to store acceleration (vector - x and y direction) at one time step
    double a[2];	   
    /// Variable to store density at one time step
    double rho;
    /// Variable to store Pressure at one time step
    double P;


    /// Variable to store previous steps position (vector - x and y direction)
    double prev_x[2];
    /// Variable to store previous steps velocity (vector - x and y direction)
    double prev_v[2];
    /// Variable to store previous density
    double prev_rho = 1000.;
    /// Variable to store D variable on solver at one time step
    double D;
    /// Variable to store Acceleration variable (vector - x and y direction) on solver at start of time step
    double gravity[2];

    /// Marker for each particles to indicate whether this particle lies on boundary or not
    bool boundary;

    /// Provides information about domain to individual particles
    static SPH_main* main_data;		

    /// Index in neighbour finding array
    int list_num[2];
    /// Index in particle_list
    int index = -1;        
    /// Index within current gridcell
    int grid_index = -1;
};

/// Class object that contain domain / system information
///
/// This domain has several methods to simulate SPH
class SPH_main
{
public:
    /// Constructor for SPH class
    SPH_main();

    /// Setting up particles in the sytem
    ///
    /// Both boundary and non boundary particles are set up in this method
    /// @return void
    /// @param min minimum value for the domain
    /// @param max maximum value for the domain
    /// @param geometry "default": default geometry; "wave": wave shoaling on beach; "step-wise": fluid flowing down steps; "bubble": water flowing over fixed spheres; "floor": fluid flowing through gap in the floor
    /// @attention Please take care when setting up the boundary. Note that boundary conditions have only been implemented for default and wave geometries.
    void place_points(double* min, double* max, string geometry);

    /// Initializing constant value and creating particle class in the domain
    ///
    /// This is the place where you set up your system length and width.
    /// @return void
    /// @see dx
    /// @see c0
    /// @see h_fac
    /// @see h
    /// @see rho0
    /// @see mu
    /// @see constant_kernel
    /// @see m0
    /// @see gama
    /// @attention Please pay attention that we are using International System of Units
    void set_values(void);

    /// Initializing constant value and creating particle class in the domain
    /// @return void
    void initialise_grid(void);

    /// Allocating all the points to the search grid (assumes that index has been appropriately updated)
    ///
    /// This function will be called to update particle grid location for each time step
    /// @return void
    void allocate_to_grid(void);

    /// Iterating neighbour cells wihtin 2h radius for every particles 
    ///
    /// This function will be called to gather information in SPH simulation 
    /// @param *part a single particle in particle list
    /// @param dt current dt
    /// @return void
    void neighbour_iterate(SPH_particle* part, double dt); 

    /// Calculating smooth kernel
    /// @param r distance between two particles (double type)
    /// @return a value (double type) which represent smoothing kernel parameter
    double W(double r);

    /// Calculating smooth kernel differential
    /// @param r distance between two particles (double type)
    /// @return a value (double type) which represent smoothing kernel parameter
    double differential_W(double r);

    /// Smoothing density distribution in the field to avoid extreme discrepancy between particles
    /// @param *part a single particle in particle list
    /// @param dt current dt
    /// @return void
    void smooth(SPH_particle* part, double dt);

    /// Updating particles variable using forward euler function
    /// @param *part a single particle in particle list
    /// @param dt current dt
    /// @return void
    void forward_euler(SPH_particle* part, double dt, string geometry);

    /// Performing half timestep using predictor corrector function
    /// @param particle_list A list of every single particle
    /// @param dt time step interval 
    /// @return void
    void predictor_corrector_half(SPH_particle* part, double dt);

    /// Updating particle variables using predictor corrector function
    /// @param particle_list A list of every single particle
    /// @param dt time step interval 
    /// @return void
    void predictor_corrector_full(SPH_particle* part, double dt);

    /// Calculating repulsive force inside the bounndaty
    /// @param dist length between two particles
    /// @return double the amount of force added into the particle
    double repulsive_f(double dist);

    /// Calculating next time step interval for adaptive timestepping
    /// @note Values for stability conditions are calculated in other functions.
    /// @return void
    void get_dt();

    /// Smoothing length
    double h;				
    /// Smoothing factor
    double h_fac;
    /// Intial Particle spacing
    double dx;

    //double time[3];

    /// Viscosity (Pa. S)
    double mu;
    /// Initial Desntiy (Kg / m^3)
    double rho0;
    /// Initial Mass (Kg)
    double m0;
    //bool boundary;
    /// Constant kernel number for W and dW/dr
    double constant_kernel;

    /// Variable to compute required dt for CFL stability condition
    double dt_CFL;
    /// Variable to compute required dt for force stability condition
    double dt_F;
    /// Variable to compute required dt for density stability condition
    double dt_A;
    /// time interval between time step
    double dt;

    /// Speed of sound (m / s)
    double c0;
    /// Gamma factor
    double gamma;

    // Dimensions of simulation region
    double min_x[2], max_x[2];				
    // Grid dimensions
    int max_list[2];

    /// List of all the particles
    vector<SPH_particle> particle_list;						
    /// Address for each particle in grid. Outer 2 are the grid, inner vector is the list of pointers in each cell
    vector<vector<vector<SPH_particle*> > > search_grid;		
};

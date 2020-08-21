#include "../includes/SPH_2D.h"
#include <cmath>
#include <complex>
#include <math.h>
#include <omp.h>
#include <iomanip>
#define _USE_MATH_DEFINES

using namespace std;

SPH_main* SPH_particle::main_data;

SPH_particle::SPH_particle()
{
    // Initialize variables for all particles
    for (int i = 0; i < 2; i++) {
        x[i] = 0.0;
        v[i] = 0.0;

        prev_x[i] = 0.0;
        prev_v[i] = 0.0;
    }
    a[0] = 0.0;
    a[1] = -9.81;

    rho = main_data->rho0;
    prev_rho = main_data->rho0;
    P = 0.0;
    D = 0.0;
    gravity[0] = 0.0;
    gravity[1] = -9.81;
    boundary = false;
}

// Get grid indices
void SPH_particle::calc_index(void)
{
    for (int i = 0; i < 2; i++)
        list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

SPH_main::SPH_main()
{
    SPH_particle::main_data = this;
}

// Set simulation parameters
void SPH_main::set_values(void)
{
    // Initialise the parameters for the domain and the simulation parameters
    min_x[0] = 0.0;
    min_x[1] = 0.0;

    max_x[0] = 20.0;
    max_x[1] = 10.0;

    dx = 0.2; // resolution in x, y
    c0 = 20.; // numerical speed of sound

    h_fac = 1.3;
    h = dx * h_fac;

    SPH_particle particle;

    // Physical properties
    rho0 = 1000.001; // kg/m^3
    mu = 0.001; // Pa s
    constant_kernel = 10. / (7 * 3.1415926 * pow(h, 2));
    m0 = pow(dx, 2) * rho0;
    gamma = 7.;

    // Initial estimates for timestep
    dt = 1.3 * h / c0;
    dt_CFL = 1.3 * h / c0;
    dt_F = 1.3 * h / c0;
    dt_A = 1.3 * h / c0;
}

// Create the grid for simulation
void SPH_main::initialise_grid(void)
{
    // Initalise simulation grid parameters
    for (int i = 0; i < 2; i++)
    {
        // Add buffer for virtual wall particles
        // Inner xmin is point 0, 0
        min_x[i] -= 2.0 * h;
        // Inner xmax is point 20, 10
        max_x[i] += 2.0 * h;

        max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
    }
    search_grid.resize(max_list[0]);
    for (int i = 0; i < max_list[0]; i++)
        search_grid[i].resize(max_list[1]);
}

// Create initial condition
void SPH_main::place_points(double* min, double* max, string geometry)
{
    // Place border and fluid points in domain initialised to a damn break scenario
    // domain boarder
    double outer_min_x = min[0];
    double outer_max_x = max[0];
    double outer_min_y = min[1];
    double outer_max_y = max[1];

    // Domain for fluid
    double boarder = 2.0 * h;
    double inner_min_x = min[0] + boarder;
    double inner_max_x = max[0] - boarder;
    double inner_min_y = min[1] + boarder;
    double inner_max_y = max[1] - boarder;

    // Bottom rectangle height
    double bottom_rectangle_height = 2.;
    double rectangle_cord_max_y = inner_min_y + bottom_rectangle_height;

    // Left square width & height
    double left_square_width = 3.;
    double left_square_height = 3.;
    double square_cord_max_x = inner_min_x + left_square_width;
    double square_cord_min_y = inner_min_y + bottom_rectangle_height;
    double square_cord_max_y = inner_min_y + bottom_rectangle_height + left_square_height;
    SPH_particle particle;

    if (geometry == "default")
    {
        ///////////////////// Red Area /////////////////////
        // Add fluid particles
        // Bottom rectangle
        // X[0] - X[19]
        // Y[0] - Y[1]
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if ((i > inner_min_x&& i < inner_max_x) && (j > inner_min_y&& j < rectangle_cord_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }
        // Left square
        // X[0] - X[2]
        // Y[2] - Y[4]
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if ((i > inner_min_x&& i < square_cord_max_x) && (j > square_cord_min_y&& j < square_cord_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }

        ///////////////////// Blue Area /////////////////////
        // Add boundary particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if (!(i > inner_min_x&& i < inner_max_x) || !(j > inner_min_y&& j < inner_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }
    }

    else if (geometry == "step-wise")
    {
        ///////////////////// Red Area /////////////////////
        // Add fluid particles

        // right-top particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if ((i > (inner_min_x + 8) && i < (inner_max_x - dx)) &&
                    (j > (inner_min_y + 8) && j < (inner_max_y - dx))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }

        ///////////////////// Blue Area /////////////////////
        // Add boundary particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                // retangle area
                if (!(i > inner_min_x&& i < inner_max_x) || !(j > inner_min_y&& j < inner_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // right-side_1
                else if ((i > (inner_min_x + 5) && i < inner_max_x) && (j > (inner_min_y + 6) && j < (inner_min_y + 7))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // left-side_1
                else if ((i > (inner_min_x) && i < inner_max_x - 5) && (j > (inner_min_y + 4) && j < (inner_min_y + 5))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // right-side_2
                else if ((i > (inner_min_x + 5) && i < inner_max_x) && (j > (inner_min_y + 2) && j < (inner_min_y + 3))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // left-side_2
                else if ((i > (inner_min_x) && i < inner_max_x - 5) && (j > (inner_min_y) && j < (inner_min_y + 1))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }
    }

    else if (geometry == "bubble")
    {
        ///////////////////// Red Area /////////////////////
        // Add fluid particles

        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                // left-top particles
                if ((i > (inner_min_x) && i < (inner_max_x - 17.5)) && (j > (inner_min_y + 8 + dx) && j < (inner_max_y - dx))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // right-top particles
                else if ((i > (inner_min_x + 7.5) && i < (inner_max_x - 5)) && (j > (inner_min_y + 8 + dx) && j < (inner_max_y - dx))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // right-top particles
                else if ((i > (inner_min_x + 17.5) && i < (inner_max_x)) && (j > (inner_min_y + 8 + dx) && j < (inner_max_y - dx))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }

        ///////////////////// Blue Area /////////////////////
        // Add boundary particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {

                // top circle
                double t1 = pow((i - 4), 2) + pow((j - 6), 2);
                double t2 = pow((i - 8), 2) + pow((j - 6), 2);
                double t3 = pow((i - 12), 2) + pow((j - 6), 2);
                double t4 = pow((i - 16), 2) + pow((j - 6), 2);
                // middle circle
                double m1 = pow((i - 2), 2) + pow((j - 4), 2);
                double m2 = pow((i - 6), 2) + pow((j - 4), 2);
                double m3 = pow((i - 10), 2) + pow((j - 4), 2);
                double m4 = pow((i - 14), 2) + pow((j - 4), 2);
                double m5 = pow((i - 18), 2) + pow((j - 4), 2);
                // bottom circle
                double b1 = pow((i - 4), 2) + pow((j - 2), 2);
                double b2 = pow((i - 8), 2) + pow((j - 2), 2);
                double b3 = pow((i - 12), 2) + pow((j - 2), 2);
                double b4 = pow((i - 16), 2) + pow((j - 2), 2);

                // retangle area boundary
                if (!(i > inner_min_x&& i < inner_max_x) || !(j > inner_min_y&& j < inner_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // left-top phase
                else if ((i > (inner_min_x) && i < (inner_max_x - 17.5)) && (j > (inner_min_y + 7) && j < (inner_max_y - 2))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // right-top phase
                else if ((i > (inner_min_x + 17.5) && i < inner_max_x) && (j > (inner_min_y + 7) && j < (inner_max_y - 2))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // top circle boundary
                else if (t1 <= 1 || t2 <= 1 || t3 <= 1 || t4 <= 1) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // middle circle boundary
                else if (m1 <= 1 || m2 <= 1 || m3 <= 1 || m4 <= 1 || m5 <= 1) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                // bottom circle boundary
                else if (b1 <= 1 || b2 <= 1 || b3 <= 1 || b4 <= 1) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }
    }

    else if (geometry == "wave") {
        double interval = 0.01 * dx;
        ///////////////////// Blue Area /////////////////////
        // Add boundary particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                double bnd = 4. / 7. * i - 40. / 7.;
                if (!(i > inner_min_x&& i < inner_max_x) || !(j > inner_min_y&& j < inner_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (10 + inner_min_x) && i < (17 + inner_min_x)) && (j < bnd)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (17 + inner_min_x) && i < inner_max_x) && (j > inner_min_y&& j < (4 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }

        ///////////////////// Red Area /////////////////////
        // Add fluid particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                double bnd = 4. / 7. * i - 40. / 7.;
                if ((i > (inner_min_x + interval) && i < (10 + inner_min_x)) && (j > (inner_min_y + interval) && j < (4 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (10 + inner_min_x) && i < (17 + inner_min_x - interval)) && (j > bnd + interval) && (j < 4)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (inner_min_x + interval) && i < (inner_max_x - 3 - interval)) && (j > (4 + inner_min_y) && j < (5 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (inner_max_x - 3 - interval) && i < (inner_max_x - interval * 2)) && (j > (4 + inner_min_y + interval) &&
                    j < (5 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (1 + inner_min_x) && i < (3 + inner_min_x)) && (j < (i + 4 + inner_min_y) && j < (-i + 8 + inner_min_y) && j>5)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
            }
        }
    }

    else if (geometry == "floor")
    {
        double interval = 0.01 * dx;
        ///////////////////// Blue Area /////////////////////
        // Add boundary particles
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if (!(i > inner_min_x&& i < inner_max_x) || !(j > inner_min_y&& j < inner_max_y)) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > inner_min_x&& i < (5 + inner_min_x)) && (j > (5 + inner_min_y) && j < (6 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }
                else if ((i > (15 + inner_min_x) && i < (20 + inner_min_x)) && (j > (5 + inner_min_y) && j < (6 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = true;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }

            }
        }

        ///////////////////// Red Area /////////////////////
        for (double i = outer_min_x; i < outer_max_x; i += dx) {
            for (double j = outer_min_y; j < outer_max_y; j += dx) {
                if ((i > inner_min_x&& i < (20 + inner_min_x)) && (j > (8 + inner_min_y) && j < (9 + inner_min_y))) {
                    particle.x[0] = i;
                    particle.x[1] = j;
                    particle.boundary = false;
                    particle.calc_index();
                    particle_list.push_back(particle);
                }

            }
        }
    }
}

// Allocate particles to grid
void SPH_main::allocate_to_grid(void)
{
    for (int i = 0; i < max_list[0]; i++)
        for (int j = 0; j < max_list[1]; j++)
            search_grid[i][j].clear();

    for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
    {
        // Reset acceleration and D
        particle_list[cnt].a[0] = 0.;
        particle_list[cnt].a[1] = -9.8;
        particle_list[cnt].D = 0.;

        // Only allocate particles with valid grid indices
        if (particle_list[cnt].list_num[0] >= 0 && particle_list[cnt].list_num[0] < max_list[0])
        {
            if (particle_list[cnt].list_num[1] >= 0 && particle_list[cnt].list_num[1] < max_list[1])
                // passing particle by reference -- not a copy!
            {
                search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
                particle_list[cnt].index = cnt;
                particle_list[cnt].grid_index = search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].size() - 1;
            }
            // Delete particles with invalid indices
            else particle_list.erase(particle_list.begin() + cnt);
        }
        else particle_list.erase(particle_list.begin() + cnt);
    }
}

// Neighbor function iterating through half the neighbors and updating a, D for the pairs
void SPH_main::neighbour_iterate(SPH_particle* part, double dt)
{
    SPH_particle* other_part;
    double dist;			//distance between particles
    double dn[2];			//vector from 1st to 2nd particle

    // same grid cell
    if (part->list_num[0] >= 0 && part->list_num[0] < max_list[0])
        if (part->list_num[1] >= 0 && part->list_num[1] < max_list[1])
        {
            int ii = part->list_num[0];
            int jj = part->list_num[1];
            int kk = part->grid_index;
            // loop over particles with higher indices in same grid
            for (unsigned int cnt = kk + 1; cnt < search_grid[ii][jj].size(); cnt++)
            {
                other_part = search_grid[ii][jj][cnt];

                //Calculates the distance between potential neighbours
                for (int n = 0; n < 2; n++)
                    dn[n] = part->x[n] - other_part->x[n];
                dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

                if (dist < 2. * h)					//only particle within 2h
                {
                    double diff_W = differential_W(dist);
                    double v_ij[2];

                    // Acceleration
                    for (int n = 0; n < 2; n++)
                    {
                        v_ij[n] = part->v[n] - other_part->v[n];

                        part->a[n] += -m0 * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * diff_W * dn[n] / dist + \
                            mu * m0 * (1.0 / (part->rho * part->rho) + 1.0 / (other_part->rho * other_part->rho)) * diff_W * \
                            v_ij[n] / dist;

                        other_part->a[n] += -m0 * (other_part->P / (other_part->rho * other_part->rho) + part->P / (part->rho * part->rho)) * diff_W * (-dn[n]) / dist + \
                            mu * m0 * (1.0 / (other_part->rho * other_part->rho) + 1.0 / (part->rho * part->rho)) * diff_W * \
                            (-v_ij[n]) / dist;
                    }
                    // D
                    part->D += m0 * diff_W * (v_ij[0] * dn[0] + v_ij[1] * dn[1]) / dist;
                    other_part->D += m0 * diff_W * ((-v_ij[0]) * (-dn[0]) + (-v_ij[1]) * (-dn[1])) / dist;

                    // Finding dt satisfying CFL stability criterion
                    double ij_dt_CFL = h / sqrt(v_ij[0] * v_ij[0] + v_ij[1] * v_ij[1]);
                    if (ij_dt_CFL < dt_CFL)
                        dt_CFL = ij_dt_CFL;
                }
            }
        }

    // middle right cell
    if (part->list_num[0] + 1 >= 0 && part->list_num[0] + 1 < max_list[0])
        if (part->list_num[1] >= 0 && part->list_num[1] < max_list[1])
        {
            int i = part->list_num[0] + 1;
            int j = part->list_num[1];

            for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
            {
                other_part = search_grid[i][j][cnt];

                if (part != other_part)					//stops particle interacting with itself
                {
                    //Calculates the distance between potential neighbours
                    for (int n = 0; n < 2; n++)
                        dn[n] = part->x[n] - other_part->x[n];
                    dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

                    if (dist < 2. * h)					//only particle within 2h
                    {
                        double diff_W = differential_W(dist);
                        double v_ij[2];

                        // acceleration
                        for (int n = 0; n < 2; n++)
                        {
                            v_ij[n] = part->v[n] - other_part->v[n];

                            part->a[n] += -m0 * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * diff_W * dn[n] / dist + \
                                mu * m0 * (1.0 / (part->rho * part->rho) + 1.0 / (other_part->rho * other_part->rho)) * diff_W * \
                                v_ij[n] / dist;

                            other_part->a[n] += -m0 * (other_part->P / (other_part->rho * other_part->rho) + part->P / (part->rho * part->rho)) * diff_W * (-dn[n]) / dist + \
                                mu * m0 * (1.0 / (other_part->rho * other_part->rho) + 1.0 / (part->rho * part->rho)) * diff_W * \
                                (-v_ij[n]) / dist;
                        }
                        // D
                        part->D += m0 * diff_W * (v_ij[0] * dn[0] + v_ij[1] * dn[1]) / dist;
                        other_part->D += m0 * diff_W * ((-v_ij[0]) * (-dn[0]) + (-v_ij[1]) * (-dn[1])) / dist;

                        // Finding dt satisfying CFL stability criterion
                        double ij_dt_CFL = h / sqrt(v_ij[0] * v_ij[0] + v_ij[1] * v_ij[1]);
                        if (ij_dt_CFL < dt_CFL)
                            dt_CFL = ij_dt_CFL;
                    }
                }
            }
        }

    // top row
    for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
        if (i >= 0 && i < max_list[0])
        {
            int j = part->list_num[1] + 1;
            if (j >= 0 && j < max_list[1])
                for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                {
                    other_part = search_grid[i][j][cnt];

                    if (part != other_part)					//stops particle interacting with itself
                    {
                        //Calculates the distance between potential neighbours
                        for (int n = 0; n < 2; n++)
                            dn[n] = part->x[n] - other_part->x[n];
                        dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

                        if (dist < 2. * h)					//only particle within 2h
                        {
                            double diff_W = differential_W(dist);
                            double v_ij[2];

                            // acceleration
                            for (int n = 0; n < 2; n++)
                            {
                                v_ij[n] = part->v[n] - other_part->v[n];

                                part->a[n] += -m0 * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * diff_W * dn[n] / dist + \
                                    mu * m0 * (1.0 / (part->rho * part->rho) + 1.0 / (other_part->rho * other_part->rho)) * diff_W * \
                                    v_ij[n] / dist;

                                other_part->a[n] += -m0 * (other_part->P / (other_part->rho * other_part->rho) + part->P / (part->rho * part->rho)) * diff_W * (-dn[n]) / dist + \
                                    mu * m0 * (1.0 / (other_part->rho * other_part->rho) + 1.0 / (part->rho * part->rho)) * diff_W * \
                                    (-v_ij[n]) / dist;
                            }
                            // D
                            part->D += m0 * diff_W * (v_ij[0] * dn[0] + v_ij[1] * dn[1]) / dist;
                            other_part->D += m0 * diff_W * ((-v_ij[0]) * (-dn[0]) + (-v_ij[1]) * (-dn[1])) / dist;

                            // Finding dt satisfying CFL stability criterion
                            double ij_dt_CFL = h / sqrt(v_ij[0] * v_ij[0] + v_ij[1] * v_ij[1]);
                            if (ij_dt_CFL < dt_CFL)
                                dt_CFL = ij_dt_CFL;
                        }
                    }
                }
        }
}

// Get next timestep
void SPH_main::get_dt()
{
    // CFL value - adjust for stability vs speed (0.1 - 0.3)
    double c_CFL = 0.15;

    // Find minimum required dt for stability
    if (this->dt_F < this->dt_CFL)
        this->dt = c_CFL * this->dt_F;
    else
    {
        if (this->dt_A < this->dt_F)
            this->dt = c_CFL * this->dt_A;
        else this->dt = c_CFL * this->dt_CFL;
    }

    //// PRINT EVERY TIMESTEP INTERVAL (for debugging) ////
    //cout << setprecision(5) << fixed;
    //cout << this->dt << endl;

    // Reset to initial underestimates for timestep
    this->dt_CFL = h / c0;
    this->dt_F = h / c0;
    this->dt_A = h / c0;
}

double SPH_main::W(double r) 
{
    double w = 0.0;
    double q = abs(r) / h;

    if (q < 1)   w = 1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3);
    else if (q >= 1 && q <= 2)   w = 0.25 * pow(2 - q, 3);
    else    w = 0;
    return w * constant_kernel;
}

double SPH_main::differential_W(double r) 
{
    double diff_w = 0.0;
    double q = abs(r) / h;

    if (q <= 1)   diff_w = -3. * q + 9.0 / 4.0 *q * q;
    else    diff_w = -0.75 * (2 - q) * (2 -q);
    return diff_w * (constant_kernel / h);
}

// Density smoothing function
void SPH_main::smooth(SPH_particle* part, double dt)
{
    SPH_particle* other_part;
    double dist;			//distance between particles
    double dn[2];			//vector from 1st to 2nd particle

    double up = 0;
    double down = 0;

    for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
        if (i >= 0 && i < max_list[0])
            for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
                if (j >= 0 && j < max_list[1])
                {
                    for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                    {
                        other_part = search_grid[i][j][cnt];

                        // Calculates the distance between potential neighbours
                        for (int n = 0; n < 2; n++)
                            dn[n] = part->x[n] - other_part->x[n];

                        dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                        double temp = W(dist);

                        // Only consider particles within 2h distance
                        if (dist <= 2. * h)
                        {
                            up += temp;
                            down += temp / other_part->prev_rho;
                        }
                    }
                }
        part->rho = up / down;
}

// Repulsive boundaries
double SPH_main::repulsive_f(double dist)
{
    return pow(0.5 * this->dx / dist, 4) * (1 - dist / (0.5 * this->dx));
}

// Predictor-Corrector scheme implementation
void SPH_main::predictor_corrector_half(SPH_particle* part, double dt) {
    // Half-step
    if (!part->boundary) {
        for (int i = 0; i < 2; i++) 
        {
            part->prev_x[i] = part->x[i];
            part->prev_v[i] = part->v[i];

            // Repulsive boundaries 
            if (part->x[i] < 0.4 * dx)
                part->a[i] += (0.05 * 9.81 * 5 / dx * repulsive_f(part->x[i])) * 0.25;
            if (part->x[i] > max_x[i] - 2 * h - 0.4 * dx)
                part->a[i] -= (0.05 * 9.81 * 5 / dx * repulsive_f(max_x[i] - 2 * h - part->x[i])) * 0.25;
            
            // Update x and v for fluid particles
            part->x[i] = part->x[i] + 0.5 * dt * part->v[i];
            part->v[i] = part->v[i] + 0.5 * dt * part->a[i];
        }
        part->calc_index();
    }
    // Update density for all particles
    part->prev_rho = part->rho;
    part->rho = part->rho + 0.5 * dt * part->D;
    part->P = rho0 * pow(c0, 2) / gamma * (pow((part->rho / rho0), gamma) - 1);

        // Finding dt_F needed to satisfy force condition for stability
    double i_dt_F = sqrt(h / sqrt(part->a[0] * part->a[0] + part->a[1] * part->a[1]));
    if (i_dt_F < this->dt_F)
        this->dt_F = i_dt_F;

    // Finding dt_A needed to satisfy density condition for stability
    double i_dt_A = h / (c0 * pow(part->rho / rho0, (gamma - 1) / 2));
    if (i_dt_A < this->dt_A)
        this->dt_A = i_dt_A;
}

void SPH_main::predictor_corrector_full(SPH_particle* part, double dt) {
    // Full-step
    if (!part->boundary) {
        for (int i = 0; i < 2; i++) 
        {
            // Repulsive boundaries
            if (part->x[i] < 0.4 * dx)
                part->a[i] += (0.05 * 9.81 * 5 / dx * repulsive_f(part->x[i])) * 0.25;
            if (part->x[i] > max_x[i] - 2 * h - 0.4 * dx)
                part->a[i] -= (0.05 * 9.81 * 5 / dx * repulsive_f(max_x[i] - 2 * h - part->x[i])) * 0.25;

            // Update x and v for fluid particles
            double x_ = part->prev_x[i] + 0.5 * dt * part->v[i];
            part->x[i] = 2 * x_ - part->prev_x[i];
            double v_ = part->prev_v[i] + 0.5 * dt * part->a[i];
            part->v[i] = 2 * v_ - part->prev_v[i];
        }
        part->calc_index();
    }
    // Update density for all particles
    double rho_ = part->prev_rho + 0.5 * dt * part->D;
    part->rho = 2 * rho_ - part->prev_rho;
    part->P = rho0 * pow(c0, 2) / gamma * (pow((part->rho / rho0), gamma) - 1);
}

// Forward euler scheme implementation
void SPH_main::forward_euler(SPH_particle* part, double dt, string geometry)
{
    if (geometry == "default")
    {
        if (part->boundary == false)
        {
            for (int n = 0; n < 2; n++)
            {
                // Boundary forces
                if (part->x[n] < 0.4 * dx)
                    part->a[n] += (0.05 * 9.81 * 5 / dx * repulsive_f(part->x[n])) * 0.25;
                if (part->x[n] > max_x[n] - 2 * h - 0.4 * dx)
                    part->a[n] -= (0.05 * 9.81 * 5 / dx * repulsive_f(max_x[n] - 2 * h - part->x[n])) * 0.25;

                // Update position and velocity
                part->x[n] += dt * part->v[n];
                part->v[n] += dt * part->a[n];
            }
        }
    }

    // For beach wave initial condition
    else if (geometry == "wave")
    {
        if (part->boundary == false)
            for (int n = 0; n < 2; n++)
            {
                if (part->x[n] < 0.5 * dx)
                    part->a[n] += 9.81 * 5 / dx * repulsive_f(part->x[n]);

                if (part->x[n] > max_x[n] - 2 * h - 0.5 * dx)
                    part->a[n] -= 9.81 * 5 / dx * repulsive_f(max_x[n] - 2 * h - part->x[n]);
            }

            if (part->x[0] >= 17 && part->x[1] > 4 && part->x[1] < 4 + 0.5 * dx)
                part->a[1] += 9.81 * 5 / dx * repulsive_f(part->x[1] - 4);

            double dist_to_shoaling = abs(4. * part->x[0] - 7. * part->x[1] - 40.) / sqrt(4 * 4 + 7 * 7);
            
            if (dist_to_shoaling < 0.5 * dx && part->x[0] < 17)
                for (int i = 0; i < 2; i++)
                {
                    part->v[i] = -1. * part->v[i];
                    part->x[i] += dt * part->v[i];
                    part->prev_x[i] = part->x[i];
                }

            else
                for (int n = 0; n < 2; n++)
                {
                    part->x[n] += dt * part->v[n];
                    part->prev_x[n] = part->x[n];
                    part->v[n] += dt * part->a[n];
                }
    }

    // Update density and pressure
    part->rho += dt * part->D;
    part->P = rho0 * pow(c0, 2) / gamma * (pow((part->rho / rho0), gamma) - 1);
    part->prev_rho = part->rho;

    part->calc_index();

    // Finding dt_F needed to satisfy force condition for stability
    double i_dt_F = sqrt(h / sqrt(part->a[0] * part->a[0] + part->a[1] * part->a[1]));
    if (i_dt_F < this->dt_F)
        this->dt_F = i_dt_F;

    // Finding dt_A needed to satisfy density condition for stability
    double i_dt_A = h / (c0 * pow(part->rho / rho0, (gamma - 1) / 2));
    if (i_dt_A < this->dt_A)
        this->dt_A = i_dt_A;
}

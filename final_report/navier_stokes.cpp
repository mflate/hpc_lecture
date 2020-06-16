#include <iostream>
#include <vector>

typedef std::vector<std::vector<double>> matrix2d;

struct simulation_t {
    int nx;
    int ny;
    int nt;
    int nit;
    double dx;
    double dy;
    double dt;
};

void build_up_b(matrix2d &b, int rho, matrix2d u, matrix2d v, simulation_t sim) {
    
    double dx = sim.dx;
    double dy = sim.dy;
    double dt = sim.dt;

    // Calculate b, which is a part of the pressure field formula
    for (int j = 1; j < sim.ny - 1; j++) {
        for (int i = 1; i < sim.nx - 1; i++) {
            b[j][i] =   rho * (1 / dt * ((u[j][i+1] - u[j][i-1])/(2 * dx) + (v[j+1][i] - v[j-1][i])/(2 * dy)) -
                                   (u[j][i+1] - u[j][i-1])*(u[j][i+1] - u[j][i-1])/(2*dx)/(2*dx)    -
                                   (u[j][i+1] - u[j][i-1])/(2*dx)*(v[j+1][i] - v[j-1][i])/(2*dy)    -
                                   (v[j+1][i] - v[j-1][i])*(v[j+1][i] - v[j-1][i])/(2*dy)/(2*dy));
        }
    }
}

void pressure_poisson(matrix2d &p, matrix2d b, simulation_t sim) {

    double dx = sim.dx;
    double dy = sim.dy;

    matrix2d pn(sim.ny, std::vector<double> (sim.nx, 0));

    //  Calculate pressure field
    for (int k = 0; k < sim.nit; k++) {

            // Make copy of p
        for (int j = 0; j < sim.ny; j++) {
            for (int i = 0; i < sim.nx; i++) {
                pn[j][i] = p[j][i];
            }
        }

        for (int j = 1; j < sim.ny - 1; j++) {
            for (int i = 1; i < sim.nx - 1; i++) {
                p[j][i] = ((pn[j][i+1] + pn[j][i-1])*dy*dy + (pn[j+1][i] + pn[j-1][i])*dx*dx)/(2*(dx*dx + dy*dy)) - 
                dx*dx*dy*dy/(2 * (dx*dx + dy*dy))*b[j][i];
            }
        }

        // Boundery conditions
        for (int j = 0; j < sim.ny; j++) p[j][sim.ny-1] = p[j][sim.ny-2]; 
        for (int i = 0; i < sim.nx; i++) p[0][i]    = p[1][i];  
        for (int j = 0; j < sim.ny; j++) p[j][0]    = p[j][1];
        for (int i = 0; i < sim.nx; i++) p[sim.nx-1][i] = 0;
    }
}


void cavity_flow(matrix2d &u, matrix2d &v, matrix2d &p, int rho, double nu, simulation_t sim){

    double dx = sim.dx;
    double dy = sim.dy;
    double dt = sim.dt;

    matrix2d un(sim.ny, std::vector<double> (sim.nx, 0));
    matrix2d vn(sim.ny, std::vector<double> (sim.nx, 0)); 
    matrix2d b(sim.ny, std::vector<double> (sim.nx, 0));

    // Iteration over time
    for (int t = 0; t < sim.nt; t++) {

        // Make copy of velocity field
        for (int j = 0; j < sim.ny; j++) {
            for (int i = 0; i < sim.nx; i++) {
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }

        //  Calculate pressure field
        build_up_b(b, rho, u, v, sim);
        pressure_poisson(p, b, sim);

        for (int j = 1; j < sim.ny - 1; j++) {
            for (int i = 1; i < sim.nx - 1; i++) {
                
                // Velocity in x direction
                u[j][i] = un[j][i] -    un[j][i] * dt / dx * (un[j][i] - un[j][i-1]) - 
                                        vn[j][i] * dt / dx * (un[j][i] - un[j-1][i]) -
                                        dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1]) +
                                        nu * (dt / (dx*dx) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
                                              dt / (dy*dy) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]));
            
                // Velocity in y direction
                v[j][i] = vn[j][i] -    un[j][i] * dt / dx * (vn[j][i] - vn[j][i-1]) - 
                                        vn[j][i] * dt / dx * (vn[j][i] - vn[j-1][i]) -
                                        dt / (2 * rho * dy) * (p[j+1][i] - p[j-1][i]) +
                                        nu * (dt / (dx*dx) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
                                              dt / (dy*dy) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]));
            }
        }
        
        // Boundary conditions
        for (int j = 0; j < sim.ny; j++) {
            u[j][0]     = 0;
            u[j][sim.nx-1]  = 0;
            v[j][0]     = 0;
            v[j][sim.nx-1]  = 0;
        }

        for (int i = 0; i < sim.nx; i++) {
            u[0][i]     = 0;
            u[sim.ny-1][i]  = 1;
            v[0][i]     = 0;
            v[sim.ny-1][i]  = 0;
        }

    }

} 

int main() {

    int nx = 112;
    int ny = 11;

    simulation_t sim = {
        .nx = nx,
        .ny = ny,
        .nt = 100,
        .nit = 50,
        .dx = 2. /(nx - 1),
        .dy = 2. /(ny - 1),
        .dt = .001
    };
    
    int rho = 1;
    double nu = .1;

    // Initializing matrices    
    matrix2d u(ny, std::vector<double> (nx, 0));
    matrix2d v(ny, std::vector<double> (nx, 0));
    matrix2d p(ny, std::vector<double> (nx, 0));
    matrix2d b(ny, std::vector<double> (nx, 0)); 

    for (int i = 0; i < ny; i++) { 
        for (int j = 0; j < nx; j++) { 
            u[i][j] = 0;
            v[i][j] = 0;
            p[i][j] = 0;
            b[i][j] = 0;
        } 
    } 
    // Simulation
    cavity_flow(u, v, p, rho, nu, sim);

    // Print pressure matrix for comparison
    for (int i = 0; i < nx; i++) { 
        for (int j = 0; j < ny; j++) { 
            std::cout << p[i][j] << ' ';
        } 
        std::cout << std::endl;
    } 
}


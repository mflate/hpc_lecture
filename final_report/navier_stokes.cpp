#include <cstdio>
#include <vector>
#include <chrono>
#include <assert.h>


#define CHECK_ERROR

typedef std::vector<std::vector<float>> matrix2d;

struct simulation_t {
    int nx;
    int ny;
    int nt;
    int nit;
    float dx;
    float dy;
    float dt;
};

void build_up_b(matrix2d &b, int rho, matrix2d u, matrix2d v, simulation_t sim) {
    
    float dx = sim.dx;
    float dy = sim.dy;
    float dt = sim.dt;

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

    float dx = sim.dx;
    float dy = sim.dy;

    matrix2d pn(sim.ny, std::vector<float> (sim.nx, 0));

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


void cavity_flow(matrix2d &u, matrix2d &v, matrix2d &p, int rho, float nu, simulation_t sim){
    
    float dx = sim.dx;
    float dy = sim.dy;
    float dt = sim.dt;

    matrix2d un(sim.ny, std::vector<float> (sim.nx, 0));
    matrix2d vn(sim.ny, std::vector<float> (sim.nx, 0)); 
    matrix2d b(sim.ny, std::vector<float> (sim.nx, 0));

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

    int nx = 11;
    int ny = 11;

    simulation_t sim = {
        .nx = nx,
        .ny = ny,
        .nt = 100,
        .nit = 50,
        .dx = (float) 2.0 /(nx - 1),
        .dy = (float) 2.0 /(ny - 1),
        .dt = .001
    };
    
    int rho = 1;
    float nu = .1;

    // Initializing matrices    
    matrix2d u(ny, std::vector<float> (nx, 0));
    matrix2d v(ny, std::vector<float> (nx, 0));
    matrix2d p(ny, std::vector<float> (nx, 0));
    matrix2d b(ny, std::vector<float> (nx, 0)); 

    for (int j = 0; j < ny; j ++){
        for (int i = 0; i < nx; i++) {
            u[j][i] = 0;
            v[j][i] = 0;
            p[j][i] = 0;
            b[j][i] = 0;
        }
    }

    // Simulation
    auto tic = std::chrono::steady_clock::now();
    cavity_flow(u, v, p, rho, nu, sim);
    auto toc = std::chrono::steady_clock::now();

    double time = std::chrono::duration<double>(toc - tic).count();
    printf("Time: %lf s \n",time);
    
    #ifdef CHECK_ERROR
        
        assert(sim.ny == 11 && sim.nx == 11);

        float p_ref[sim.ny][sim.nx] = {
            {-7.3091e-02, -7.3091e-02, -6.0065e-02, -4.5810e-02, -2.2408e-02,  9.8650e-04, 2.4225e-02,  4.7980e-02,  6.2242e-02,  7.5588e-02,  7.5588e-02},
            {-7.3091e-02, -7.3091e-02, -6.0065e-02, -4.5810e-02, -2.2408e-02,  9.8650e-04, 2.4225e-02,  4.7980e-02,  6.2242e-02,  7.5588e-02,  7.5588e-02},
            {-4.9689e-02, -4.9689e-02, -4.0887e-02, -3.1324e-02, -1.5214e-02,  8.4575e-04, 1.6725e-02,  3.3153e-02,  4.2672e-02,  5.1742e-02,  5.1742e-02},
            {-1.1074e-01, -1.1074e-01, -8.8839e-02, -6.8636e-02, -3.2909e-02,  8.7371e-04, 3.4370e-02,  7.0752e-02,  9.0949e-02,  1.1355e-01,  1.1355e-01},
            {-4.4704e-02, -4.4704e-02, -3.5668e-02, -2.8150e-02, -1.3219e-02,  5.7227e-04, 1.4047e-02,  2.9508e-02,  3.6926e-02,  4.6455e-02,  4.6455e-02},
            {-2.0766e-01, -2.0766e-01, -1.5675e-01, -1.2354e-01, -5.6565e-02,  4.7075e-04, 5.6858e-02,  1.2526e-01,  1.5837e-01,  2.1128e-01,  2.1128e-01},
            {-4.9574e-02, -4.9574e-02, -3.6233e-02, -3.0298e-02, -1.3279e-02,  1.3012e-04, 1.2944e-02,  3.0973e-02,  3.6697e-02,  5.1263e-02,  5.1263e-02},
            {-4.2955e-01, -4.2955e-01, -2.8484e-01, -2.2962e-01, -9.7103e-02, -6.5522e-04, 9.4288e-02,  2.2939e-01,  2.8381e-01,  4.3393e-01,  4.3393e-01},
            {-5.4136e-02, -5.4136e-02, -3.1682e-02, -2.8754e-02, -1.1243e-02, -3.6575e-04, 9.4863e-03,  2.8524e-02,  3.0682e-02,  5.6854e-02,  5.6854e-02},
            {-9.6032e-01, -9.6032e-01, -4.7186e-01, -3.8849e-01, -1.4617e-01, -2.6026e-03, 1.3788e-01,  3.8276e-01,  4.6265e-01,  9.5444e-01,  9.5444e-01},
            { 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00}
        };

        
        float err = 0;
        for (int j = 0; j < ny; j ++){
            for (int i = 0; i < nx; i++) {
                err += (p_ref[j][i] - p[j][i]) * (p_ref[j][i] - p[j][i])/(sim.nx*sim.ny);
            }
        }
        printf("Error: %f\n", err);
    #endif

}

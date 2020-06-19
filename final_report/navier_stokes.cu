#include <iostream>
#include <chrono>
#include <assert.h>

using namespace std;

// Option for comparing output pressure to given python scipts
// Only works for nx x ny = 11x11 
#define CHECK_ERROR                

// Struct for storing simulation variables
struct simulation_t {
    int nx;
    int ny;
    int nt;
    int nit;
    float dx;
    float dy;
    float dt;
};


__global__ void build_up_b(
    float* b, 
    const float *u, 
    const float *v, 
    int rho,
    simulation_t sim
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Serialized index
    if (i >= sim.nx*sim.ny) return;
    
    int ip1 = i + 1;                // serialized index i + 1
    int im1 = i - 1;                // serialized index i - 1
    int jp1 = i + sim.nx;           // serialized index j + 1
    int jm1 = i - sim.nx;           // serialized index j - 1
    
    // 2d subscripts
    int i_mat = i % sim.nx;
    int j_mat = i / sim.nx;

    // Border conditions
    if (i_mat == 0 || i_mat == sim.nx - 1 || j_mat == 0 || j_mat == sim.ny - 1) {b[i] = 0; return;}

    float dx = sim.dx;
    float dy = sim.dy;
    float dt = sim.dt;
    
    // Build up b
    b[i] =  rho * (1 / dt * ((u[ip1] - u[im1])/(2 * dx) + (v[jp1] - v[jm1])/(2 * dy)) -
                            (u[ip1] - u[im1])*(u[ip1] - u[im1])/(2*dx)/(2*dx)    -
                            (u[ip1] - u[im1])/(2*dx)*(v[jp1] - v[jm1])/(2*dy)    -
                            (v[jp1] - v[jm1])*(v[jp1] - v[jm1])/(2*dy)/(2*dy));

}


__global__ void pressure_poisson(
    float * p,
    float * pn,
    const float * b,
    simulation_t sim
) {
    
    // Serialized index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i > sim.nx*sim.ny) return;

    int ip1 = i + 1;                // serialized index i + 1
    int im1 = i - 1;                // serialized index i - 1
    int jp1 = i + sim.nx;           // serialized index j + 1
    int jm1 = i - sim.nx;           // serialized index j - 1
    
    // 2d subscripts 
    int i_mat = i % sim.nx;
    int j_mat = i / sim.nx;

    float dx = sim.dx;
    float dy = sim.dy;
    
    for (int k = 0; k < sim.nit; k++) {
        
        // Copy last interation of pressure matrix
        pn[i] = p[i];
        __syncthreads();

        if (i_mat != 0 && i_mat != sim.nx - 1 && j_mat != 0 && j_mat != sim.ny - 1) {
            
            // Update pressure equation
            p[i] = ((pn[ip1] + pn[im1])*dy*dy + (pn[jp1] + pn[jm1])*dx*dx)/(2*(dx*dx + dy*dy)) - 
                dx*dx*dy*dy/(2 * (dx*dx + dy*dy))*b[i];
        
            }

        // Border conditions
        __syncthreads();
        if (j_mat == sim.ny-1)  p[i] = 0;
        if (j_mat == 0)         p[i] = p[jp1];
        __syncthreads();
        if (i_mat == sim.nx-1)  p[i] = p[im1];
        if (i_mat == 0)         p[i] = p[ip1];
        __syncthreads();
    }
   
}

__global__ void velocity_field( 
    float * u,
    float * v,
    float * un,
    float * vn,
    const float * p,
    int rho,
    float nu,
    simulation_t sim
) {
    
    // Serialized index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i > sim.nx*sim.ny) return;


    int ip1 = i + 1;                // serialized index i + 1
    int im1 = i - 1;                // serialized index i - 1
    int jp1 = i + sim.nx;           // serialized index j + 1
    int jm1 = i - sim.nx;           // serialized index j - 1
    
    // 2d subscripts
    int i_mat = i % sim.nx;
    int j_mat = i / sim.nx;

    float dx = sim.dx;
    float dy = sim.dy;
    float dt = sim.dt;

    un[i] = u[i]; vn[i] = v[i];
    __syncthreads();

    if (i_mat != 0 && i_mat != sim.nx - 1 && j_mat != 0 && j_mat != sim.ny - 1) {

        // Velocity in x direction
        u[i] = un[i] -  un[i] * dt / dx * (un[i] - un[im1]) - 
                        vn[i] * dt / dx * (un[i] - un[jm1]) -
                        dt / (2 * rho * dx) * (p[ip1] - p[im1]) +
                        nu * (dt / (dx*dx) * (un[ip1] - 2 * un[i] + un[im1]) +
                        dt / (dy*dy) * (un[jp1] - 2 * un[i] + un[jm1]));

        // Velocity in y direction
        v[i] = vn[i] -  un[i] * dt / dx * (vn[i] - vn[im1]) - 
                        vn[i] * dt / dx * (vn[i] - vn[jm1]) -
                        dt / (2 * rho * dy) * (p[jp1] - p[jm1]) +
                        nu * (dt / (dx*dx) * (vn[ip1] - 2 * vn[i] + vn[im1]) +
                        dt / (dy*dy) * (vn[jp1] - 2 * vn[i] + vn[jm1]));
    }

    // Border conditions
    __syncthreads();
    if (i_mat == 0)        { u[i] = 0; v[i] = 0; }
    if (i_mat == sim.nx-1) { u[i] = 0; v[i] = 0; }
    __syncthreads();
    if (j_mat == 0)        { u[i] = 0; v[i] = 0; }
    if (j_mat == sim.ny-1) { u[i] = 1; v[i] = 0; }
    __syncthreads();

}


void cavity_flow(
    float * u,
    float * v,
    float * p,
    int rho,
    float nu,
    simulation_t sim
) {

    // Blocks and threads
    int M = 64;
    int N = sim.nx*sim.ny;

    auto tic = chrono::steady_clock::now();

    // Declaring necessary matrices
    float* b;
    float* un; 
    float* vn;
    float* pn;

    // Allocating memory for necessary matrices 
    cudaMalloc((void **)&b,  sizeof(float) * sim.ny * sim.nx);
    cudaMalloc((void **)&un, sizeof(float) * sim.ny * sim.nx);
    cudaMalloc((void **)&vn, sizeof(float) * sim.ny * sim.nx);
    cudaMalloc((void **)&pn, sizeof(float) * sim.ny * sim.nx);

    // Loop over time
    for (int t = 0; t < sim.nt; t++) {

        // Simulation of navier stokes equation for one timestep
        build_up_b<<<(N + M - 1)/M,M>>>(b, u, v, rho, sim);
        pressure_poisson<<<(N + M - 1)/M,M>>>(p, pn, b, sim);
        velocity_field<<<(N + M - 1)/M,M>>>(u,v,un,vn,p,rho,nu,sim);

    }
    cudaDeviceSynchronize();

    // Free up memory
    cudaFree(b);
    cudaFree(un);
    cudaFree(vn);
    cudaFree(pn);    
}


int main() {

    // Size of vector field
    int nx = 11;
    int ny = 11;

    // Numerical attributes
    simulation_t sim = {
        .nx = nx,
        .ny = ny,
        .nt = 100,
        .nit = 50,
        .dx = (float) 2.0 /(nx - 1),
        .dy = (float) 2.0 /(ny - 1),
        .dt = .001
    };
    
    // Physical attributes
    int rho = 1;
    float nu = .1;

    // Output matrices stored on host
    // u: veclocity in x direction
    // v: veclocity in y direction
    // p: pressure
    float u[ny][nx];
    float v[ny][nx];
    float p[ny][nx];

    for (int j = 0; j < ny; j ++){
        for (int i = 0; i < nx; i++) {
            u[j][i] = 0;
            v[j][i] = 0;
            p[j][i] = 0;
        }
    }

    // Output matrices stored on gpu
    float* d_u; 
    float* d_v; 
    float* d_p; 

    auto tic2 = chrono::steady_clock::now();

    // Allocation on gpu
    // Note: first call to cudaMalloc initializes cuda, takes ~10x the time compared 
    // to everything else for the nx x ny = 11x11 case
    cudaMalloc((void **)&d_u, sizeof(float) * ny * nx);
    cudaMalloc((void **)&d_v, sizeof(float) * ny * nx);
    cudaMalloc((void **)&d_p, sizeof(float) * ny * nx);

    // Copying value of host variables to gpu varables
    cudaMemcpy(d_u, u, nx*ny * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v, nx*ny * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_p, p, nx*ny * sizeof(float), cudaMemcpyHostToDevice);

    auto tic = chrono::steady_clock::now();
    
    // Simulating navier stokes equation
    cavity_flow(d_u, d_v, d_p, rho, nu, sim); 
    auto toc = chrono::steady_clock::now();

    double time = chrono::duration<double>(toc - tic).count();
    printf("Time (excluding CUDA setup): %lf s \n",time);

    // Copying value of gpu variables to host varables
    cudaMemcpy(u, d_u, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(v, d_v, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(p, d_p, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_u);
    cudaFree(d_v);allocation
    cudaFree(d_p);
    
    auto toc2 = chrono::steady_clock::now();
    time = chrono::duration<double>(toc2 - tic2).count();
    printf("Time (including CUDA setup): %lf s \n",time);

    // Error checking compared to python code
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

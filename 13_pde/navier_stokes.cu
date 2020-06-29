#include <cstdlib>
#include <cstdio>
#include <cooperative_groups.h>

#define P(i, j) ((i) * nx + (j))

void allocate_2d(float *&a, int nx, int ny){
    cudaMallocManaged(&a, nx*ny*sizeof(float));
}

__global__  void build_up_b(float *b, float rho, float dt, float *u, float *v, float dx, float dy, int nx, int ny){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    if(i < 1 || i >= ny-1 || j < 1 || j >= nx-1){
        return;
    }
    b[P(i, j)] = (
        rho * (
            1 / dt * (
                (
                    u[P(i, j+1)] - u[P(i, j-1)]
                ) / (
                    2 * dx
                ) + (
                    v[P(i+1, j)] - v[P(i-1, j)]
                ) / (
                    2 * dy
                )
            ) - (
                (
                    u[P(i, j+1)] - u[P(i, j-1)]
                ) / (
                    2 * dx
                )
            ) * (
                (
                    u[P(i, j+1)] - u[P(i, j-1)]
                ) / (
                    2 * dx
                )
            ) - 2 * (
                (
                    u[P(i+1, j)] - u[P(i-1, j)]
                ) / (
                    2 * dy
                ) * (
                    v[P(i, j+1)] - v[P(i, j-1)]
                ) / (
                    2 * dx
                )
            ) - (
                (
                    v[P(i+1, j)] - v[P(i-1, j)]
                ) / (
                    2 * dy
                )
            ) * (
                (
                    v[P(i+1, j)] - v[P(i-1, j)]
                ) / (
                    2 * dy
                )
            )
        )
    );
}


__global__ void copy(float *a, float *b, int nx, int ny){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    if(i < 0 || i >= ny || j < 0 || j >= nx){
        return;
    }
    b[P(i, j)] = a[P(i, j)];

}


__global__ void pressure_poisson_step(float *p, float *pn, float dx, float dy, float *b, int nx, int ny){
    
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    if(i < 1 || i >= ny-1 || j < 1 || j >= nx-1){
        return;
    }
    p[P(i, j)] = (
        (
            (
                pn[P(i, j+1)] + pn[P(i, j-1)]
            ) * dy * dy + (
                pn[P(i+1, j)] + pn[P(i-1, j)]
            ) * dx * dx
        ) / (
            2 * (
                dx * dx + dy * dy
            )
        ) - dx * dx * dy * dy / (
            2 * (
                dx * dx + dy * dy
            )
        ) * b[P(i, j)]
    );
}

void pressure_poisson(float *p, float *pn, float dx, float dy, float *b, int nit, int nx, int ny){
    for(int i = 0;i < nit;i++){

        dim3 grid = dim3(2, 2);
        dim3 block = dim3(32, 32);
        copy<<<grid,block>>>(p, pn, nx, ny);
        pressure_poisson_step<<<grid,block>>>(p, pn, dx, dy, b, nx, ny);
        cudaDeviceSynchronize();
        for(int i = 0;i < ny;i++){
            p[P(i, nx-1)] = p[P(i, nx-2)];
        }
        for(int i = 0;i < nx;i++){
            p[P(0, i)] = p[P(1, i)];
        }
        for(int i = 0;i < ny;i++){
            p[P(i, 0)] = p[P(i, 1)];
        }
        for(int i = 0;i < nx;i++){
            p[P(ny-1, i)] = 0;
        }

    }
}

__global__ void update_u(float *u, float *v, float *un, float *vn, float dt, float dx, float dy, float rho, float nu, float *p, int nx, int ny){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    if(i < 1 || i >= ny-1 || j < 1 || j >= nx-1){
        return;
    }
    u[P(i, j)] = (
        un[P(i, j)] - un[P(i, j)] * dt / dx * (
            un[P(i, j)] - un[P(i, j-1)]
        ) - vn[P(i, j)] * dt / dy * (
            un[P(i, j)] - un[P(i-1, j)]
        ) - dt / (
            2 * rho * dx
        ) * (
            p[P(i, j+1)] - p[P(i, j-1)]
        ) + nu * (
            dt / (dx * dx) * (
                un[P(i, j+1)] - 2 * un[P(i, j)] + un[P(i, j-1)]
            ) + dt / (dy * dy) * (
                un[P(i+1, j)] - 2 * un[P(i, j)] + un[P(i-1, j)]
            )
        )
    );

}

__global__ void update_v(float *u, float *v, float *un, float *vn, float dt, float dx, float dy, float rho, float nu, float *p, int nx, int ny){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    if(i < 1 || i >= ny-1 || j < 1 || j >= nx-1){
        return;
    }
    v[P(i, j)] = (
        vn[P(i, j)] - un[P(i, j)] * dt / dx * (
            vn[P(i, j)] - vn[P(i, j-1)]
        ) - vn[P(i, j)] * dt / dy * (
            vn[P(i, j)] - vn[P(i-1, j)]
        ) - dt /  (
            2 * rho * dy
        ) * (
            p[P(i+1, j)] - p[P(i-1, j)]
        ) + nu * (
            dt / (dx * dx) * (
                vn[P(i, j+1)] - 2 * vn[P(i, j)] + vn[P(i, j-1)]
            ) + dt / (dy * dy) * (
                vn[P(i+1, j)] - 2 * vn[P(i, j)] + vn[P(i-1, j)]
            )
        )
    );

}

void boarder(float *u, float *v, int nx, int ny){
    for(int i = 0;i < nx;i++){
        u[P(0, i)] = 0;
    }
    for(int i = 0;i < ny;i++){
        u[P(i, 0)] = 0;
    }
    for(int i = 0;i < ny;i++){
        u[P(i, nx-1)] = 0;
    }
    for(int i = 0;i < nx;i++){
        u[P(ny-1, i)] = 1;
    }
    for(int i = 0;i < nx;i++){
        v[P(0, i)] = 0;
    }
    for(int i = 0;i < nx;i++){
        v[P(ny-1, i)] = 0;
    }
    for(int i = 0;i < ny;i++){
        v[P(i, 0)] = 0;
    }
    for(int i = 0;i < ny;i++){
        v[P(i, nx-1)] = 0;
    }
}

void print2d(float *a, int nx, int ny){
    for(int i = 0;i < ny;i++){
        for(int j = 0;j < nx;j++){
            printf("%g\n", a[P(i, j)]);
        }
    }
}

void cavity_flow_step(float *u, float *un, float *v, float *vn, float dt, float dx, float dy, float *p, float *pn, float *b, float rho, float nu, int nx, int ny, int nit){
    dim3 grid = dim3(2, 2);
    dim3 block = dim3(32, 32);
    build_up_b<<<grid,block>>>(b, rho, dt, u, v, dx, dy, nx, ny);
    cudaDeviceSynchronize();
    pressure_poisson(p, pn, dx, dy, b, nit, nx, ny);
    cudaDeviceSynchronize();
    copy<<<grid,block>>>(u, un, nx, ny);
    copy<<<grid,block>>>(v, vn, nx, ny);
    update_u<<<grid,block>>>(u, v, un, vn, dt, dx, dy, rho, nu, p, nx, ny);
    update_v<<<grid,block>>>(u, v, un, vn, dt, dx, dy, rho, nu, p, nx, ny);
    cudaDeviceSynchronize();
    boarder(u, v, nx, ny);
    cudaDeviceSynchronize();
}

void cavity_flow(int nt, float *u, float *un, float *v, float *vn, float dt, float dx, float dy, float *p, float *pn, float *b, float rho, float nu, int nx, int ny, int nit){
    for(int i = 0;i < nt;i++){
        cavity_flow_step(u, un, v, vn, dt, dx, dy, p, pn, b, rho, nu, nx, ny, nit);
    }
}


int main(){
    int nx = 41;
    int ny = 41;
    int nt = 700;
    int nit = 50;
    float cx = 2;
    float cy = 2;
    float dx = cx / (nx - 1);
    float dy = cy / (ny - 1);

    float rho = 1;
    float nu = 0.1;
    float dt = 0.001;
    float *u;
    allocate_2d(u, nx, ny);
    float *v;
    allocate_2d(v, nx, ny);
    float *p;
    allocate_2d(p, nx, ny);
    float *b;
    allocate_2d(b, nx, ny);
    float *un;
    allocate_2d(un, nx, ny);
    float *vn;
    allocate_2d(vn, nx, ny);
    float *pn;
    allocate_2d(pn, nx, ny);

    for(int i = 0;i < ny;i++){
        for(int j = 0;j < nx;j++){
            u[P(i, j)] = v[P(i, j)] = p[P(i, j)] = b[P(i, j)] = 0;
        }
    }
    cavity_flow(nt, u, un, v, vn, dt, dx, dy, p, pn, b, rho, nu, nx, ny, nit);
    cudaDeviceSynchronize();
    print2d(u, nx, ny);
    print2d(v, nx, ny);
    print2d(p, nx, ny);
}

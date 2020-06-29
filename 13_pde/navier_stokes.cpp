#include <cstdlib>
#include <cstdio>
#include <iostream>

void allocate_2d(float **&a, int nx, int ny){
    a = (float**)calloc(ny, sizeof(float*));
    for(int i = 0;i < ny;i++){
        a[i] = (float*)calloc(nx, sizeof(float));
    }
}

void build_up_b(float **b, float rho, float dt, float **u, float **v, float dx, float dy, int nx, int ny){
    for(int i = 1;i < ny-1;i++){
        for(int j = 1;j < nx-1;j++){
            b[i][j] = (
                rho * (
                    1 / dt * (
                        (
                            u[i][j+1] - u[i][j-1]
                        ) / (
                            2 * dx
                        ) + (
                            v[i+1][j] - v[i-1][j]
                        ) / (
                            2 * dy
                        )
                    ) - (
                        (
                            u[i][j+1] - u[i][j-1]
                        ) / (
                            2 * dx
                        )
                    ) * (
                        (
                            u[i][j+1] - u[i][j-1]
                        ) / (
                            2 * dx
                        )
                    ) - 2 * (
                        (
                            u[i+1][j] - u[i-1][j]
                        ) / (
                            2 * dy
                        ) * (
                            v[i][j+1] - v[i][j-1]
                        ) / (
                            2 * dx
                        )
                    ) - (
                        (
                            v[i+1][j] - v[i-1][j]
                        ) / (
                            2 * dy
                        )
                    ) * (
                        (
                            v[i+1][j] - v[i-1][j]
                        ) / (
                            2 * dy
                        )
                    )
                )
            );
        }
    }
}

void copy(float **a, float **b, int nx, int ny){
    for(int i = 0;i < ny;i++){
        for(int j = 0;j < nx;j++){
            b[i][j] = a[i][j];
        }
    }
}

void pressure_poisson(float **p, float **pn, float dx, float dy, float **b, int nit, int nx, int ny){
    for(int i = 0;i < nit;i++){
        copy(p, pn, nx, ny);
        for(int i = 1;i < ny-1;i++){
            for(int j = 1;j < nx-1;j++){
                p[i][j] = (
                    (
                        (
                            pn[i][j+1] + pn[i][j-1]
                        ) * dy * dy + (
                            pn[i+1][j] + pn[i-1][j]
                        ) * dx * dx
                    ) / (
                        2 * (
                            dx * dx + dy * dy
                        )
                    ) - dx * dx * dy * dy / (
                        2 * (
                            dx * dx + dy * dy
                        )
                    ) * b[i][j]
                );
            }
        }
        for(int i = 0;i < ny;i++){
            p[i][nx-1] = p[i][nx-2];
        }
        for(int i = 0;i < nx;i++){
            p[0][i] = p[1][i];
        }
        for(int i = 0;i < ny;i++){
            p[i][0] = p[i][1];
        }
        for(int i = 0;i < nx;i++){
            p[ny-1][i] = 0;
        }

    }
}

void update_u(float **u, float **v, float **un, float **vn, float dt, float dx, float dy, float rho, float nu, float **p, int nx, int ny){
    for(int i = 1;i < ny-1;i++){
        for(int j = 1;j < nx-1;j++){
            u[i][j] = (
                un[i][j] - un[i][j] * dt / dx * (
                    un[i][j] - un[i][j-1]
                ) - vn[i][j] * dt / dy * (
                    un[i][j] - un[i-1][j]
                ) - dt / (
                    2 * rho * dx
                ) * (
                    p[i][j+1] - p[i][j-1]
                ) + nu * (
                    dt / (dx * dx) * (
                        un[i][j+1] - 2 * un[i][j] + un[i][j-1]
                    ) + dt / (dy * dy) * (
                        un[i+1][j] - 2 * un[i][j] + un[i-1][j]
                    )
                )
            );
        }
    }
}

void update_v(float **u, float **v, float **un, float **vn, float dt, float dx, float dy, float rho, float nu, float **p, int nx, int ny){
    for(int i = 1;i < ny-1;i++){
        for(int j = 1;j < nx-1;j++){
            v[i][j] = (
                vn[i][j] - un[i][j] * dt / dx * (
                    vn[i][j] - vn[i][j-1]
                ) - vn[i][j] * dt / dy * (
                    vn[i][j] - vn[i-1][j]
                ) - dt /  (
                    2 * rho * dy
                ) * (
                    p[i+1][j] - p[i-1][j]
                ) + nu * (
                    dt / (dx * dx) * (
                        vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]
                    ) + dt / (dy * dy) * (
                        vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j]
                    )
                )
            );
        }
    }
}

void boarder(float **u, float **v, int nx, int ny){
    for(int i = 0;i < nx;i++){
        u[0][i] = 0;
    }
    for(int i = 0;i < ny;i++){
        u[i][0] = 0;
    }
    for(int i = 0;i < ny;i++){
        u[i][nx-1] = 0;
    }
    for(int i = 0;i < nx;i++){
        u[ny-1][i] = 1;
    }
    for(int i = 0;i < nx;i++){
        v[0][i] = 0;
    }
    for(int i = 0;i < nx;i++){
        v[ny-1][i] = 0;
    }
    for(int i = 0;i < ny;i++){
        v[i][0] = 0;
    }
    for(int i = 0;i < ny;i++){
        v[i][nx-1] = 0;
    }
}

void print2d(float **a, int nx, int ny){
    for(int i = 0;i < ny;i++){
        for(int j = 0;j < nx;j++){
            printf("%g\n", a[i][j]);
        }
    }
}

void cavity_flow_step(float **u, float **un, float **v, float **vn, float dt, float dx, float dy, float **p, float **pn, float **b, float rho, float nu, int nx, int ny, int nit){
    build_up_b(b, rho, dt, u, v, dx, dy, nx, ny);
    pressure_poisson(p, pn, dx, dy, b, nit, nx, ny);
    copy(u, un, nx, ny);
    copy(v, vn, nx, ny);
    update_u(u, v, un, vn, dt, dx, dy, rho, nu, p, nx, ny);
    update_v(u, v, un, vn, dt, dx, dy, rho, nu, p, nx, ny);
    boarder(u, v, nx, ny);
}

void cavity_flow(int nt, float **u, float **un, float **v, float **vn, float dt, float dx, float dy, float **p, float **pn, float **b, float rho, float nu, int nx, int ny, int nit){
    for(int i = 0;i < nt;i++){
        cavity_flow_step(u, un, v, vn, dt, dx, dy, p, pn, b, rho, nu, nx, ny, nit);
    }
}


int main(){
    int nx = 41;
    int ny = 41;
    int nt = 3;
    int nit = 50;
    float cx = 2;
    float cy = 2;
    float dx = cx / (nx - 1);
    float dy = cy / (ny - 1);

    float rho = 1;
    float nu = 0.1;
    float dt = 0.001;
    float **u;
    allocate_2d(u, nx, ny);
    float **v;
    allocate_2d(v, nx, ny);
    float **p;
    allocate_2d(p, nx, ny);
    float **b;
    allocate_2d(b, nx, ny);
    float **un;
    allocate_2d(un, nx, ny);
    float **vn;
    allocate_2d(vn, nx, ny);
    float **pn;
    allocate_2d(pn, nx, ny);

    for(int i = 0;i < ny;i++){
        for(int j = 0;j < nx;j++){
            u[i][j] = v[i][j] = p[i][j] = b[i][j] = 0.;
        }
    }

    cavity_flow(nt, u, un, v, vn, dt, dx, dy, p, pn, b, rho, nu, nx, ny, nit);

    print2d(u, nx, ny);
    print2d(v, nx, ny);
    print2d(p, nx, ny);
}

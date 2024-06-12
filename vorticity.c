#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vorticity.h"

// Converts 2D coordinates to a 1D array index
#define IX(i,j) ((i)+(N+2)*(j))

// Macro to loop through each cell in the fluid grid
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

// Computes the vorticity of the velocity field
void compute_vorticity(int N, float *u, float *v, float *w) {
    int i, j;
    FOR_EACH_CELL
        w[IX(i, j)] = 0.5f * (v[IX(i + 1, j)] - v[IX(i - 1, j)] - u[IX(i, j + 1)] + u[IX(i, j - 1)]);
    END_FOR
}

// Computes the gradient of the vorticity magnitude
void compute_eta(int N, float *w, float *eta_x, float *eta_y) {
    int i, j;
    FOR_EACH_CELL
        eta_x[IX(i, j)] = 0.5f * (fabs(w[IX(i + 1, j)]) - fabs(w[IX(i - 1, j)]));
        eta_y[IX(i, j)] = 0.5f * (fabs(w[IX(i, j + 1)]) - fabs(w[IX(i, j - 1)]));
    END_FOR
}

// Normalizes the vorticity location vector
void normalize_eta(int N, float *eta_x, float *eta_y, float *norm_eta_x, float *norm_eta_y) {
    int i, j;
    float mag;
    FOR_EACH_CELL
        mag = sqrt(eta_x[IX(i, j)] * eta_x[IX(i, j)] + eta_y[IX(i, j)] * eta_y[IX(i, j)]) + 1e-5f;
        norm_eta_x[IX(i, j)] = eta_x[IX(i, j)] / mag;
        norm_eta_y[IX(i, j)] = eta_y[IX(i, j)] / mag;
    END_FOR
}

// Adds the vorticity confinement force to the velocity field
void vorticity_confinement(int N, float *u, float *v, float *w, float *eta_x, float *eta_y, float *norm_eta_x, float *norm_eta_y, float dt, float epsilon) {
    int i, j;
    float force_x, force_y;
    FOR_EACH_CELL
        force_x = epsilon * (norm_eta_y[IX(i, j)] * w[IX(i, j)]);
        force_y = epsilon * (-norm_eta_x[IX(i, j)] * w[IX(i, j)]);
        u[IX(i, j)] += dt * force_x;
        v[IX(i, j)] += dt * force_y;

        // Debug outputs
        // printf("Cell (%d, %d): Force X: %f, Force Y: %f\n", i, j, force_x, force_y);
    END_FOR
}

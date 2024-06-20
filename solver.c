#include <stdlib.h>
#include "vorticity.h"
#include <math.h>

// Converts 2D coordinates to a 1D array index
#define IX(i,j) ((i)+(N+2)*(j))

// Swaps the pointers of two float arrays
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

// Macro to loop through each cell in the fluid grid
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

// Adds a source term to the field x. 
// This represents external influences like user input.
void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2);
	// Iterate over each cell in the grid
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i]; // Add source term scaled by the time step dt
}

// Sets the boundary conditions. This ensures the fluid does not leak out of the container.
void set_bnd(int N, int b, float *x) {
    int i, j;
    // int hole_radius = N / 10; // Adjust the size of the hole
    // int hole_center_x = (N + 2) / 2;
    // int hole_center_y = (N + 2) / 2;

    // for (i = 1; i <= N; i++) {
    //     for (j = 1; j <= N; j++) {
    //         // Create a vertical wall in the middle except for the hole
    //         if (i == hole_center_x && abs(j - hole_center_y) > hole_radius) {
    //             x[IX(i, j)] = 0; // Wall
    //         }
    //     }
    // }

    for (i = 1; i <= N; i++) {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)]; // Left boundary
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)]; // Right boundary
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)]; // Bottom boundary
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)]; // Top boundary
    }

    // Corner boundaries
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}


// Solves the linear system using Gauss-Seidel relaxation to simulate diffusion
void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) { // Iteration count for relaxation
		FOR_EACH_CELL
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
		END_FOR
		set_bnd ( N, b, x ); // Apply boundary conditions
	}
}

// Handles the diffusion of the field x
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
{
	float a=dt*diff*N*N;
	lin_solve ( N, b, x, x0, a, 1+4*a );
}

// Handles the advection of the field d through the velocity field (u, v)
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    FOR_EACH_CELL
        x = i - dt0 * u[IX(i,j)];
        y = j - dt0 * v[IX(i,j)];
        
        if (x < 0.5f) {
            x = 0.5f;
        } 
        if (x > N + 0.5f) {
            x = N + 0.5f;
        }
        i0 = (int)x; 
        i1 = i0 + 1;
        
        if (y < 0.5f) {
            y = 0.5f;
        }
        if (y > N + 0.5f) {
            y = N + 0.5f;
        }
        j0 = (int)y; 
        j1 = j0 + 1;
        
        s1 = x - i0; 
        s0 = 1 - s1; 
        t1 = y - j0; 
        t0 = 1 - t1;
        
        d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) +
                     s1 * (t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]);
    END_FOR
    set_bnd ( N, b, d );
}

// Projects the velocity field to make it mass-conserving (divergence-free)
void project ( int N, float * u, float * v, float * p, float * div )
{
	int i, j;

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

	lin_solve ( N, 0, p, div, 1, 4 );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

// Executes a single time step for the density field
void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
	add_source ( N, x, x0, dt ); // Add density sources
	SWAP ( x0, x ); 
	diffuse ( N, 0, x, x0, diff, dt ); // Diffuse the density
	SWAP ( x0, x ); 
	advect ( N, 0, x, x0, u, v, dt ); // Advect the density
}

// Executes a single time step for the velocity field
void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt, float epsilon) {
    float *w, *eta_x, *eta_y, *norm_eta_x, *norm_eta_y;
    w = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
    eta_x = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
    eta_y = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
    norm_eta_x = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
    norm_eta_y = (float *)malloc((N + 2) * (N + 2) * sizeof(float));

    add_source(N, u, u0, dt); add_source(N, v, v0, dt); // Add velocity sources
    SWAP(u0, u); diffuse(N, 1, u, u0, visc, dt); // Diffuse the velocity
    SWAP(v0, v); diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0); // Make the velocity field divergence-free

    // Vorticity Confinement
    compute_vorticity(N, u, v, w);
    compute_eta(N, w, eta_x, eta_y);
    normalize_eta(N, eta_x, eta_y, norm_eta_x, norm_eta_y);
    vorticity_confinement(N, u, v, w, eta_x, eta_y, norm_eta_x, norm_eta_y, dt, epsilon);

    SWAP(u0, u); SWAP(v0, v);
    advect(N, 1, u, u0, u0, v0, dt); advect(N, 2, v, v0, u0, v0, dt); // Advect the velocity
    project(N, u, v, u0, v0); // Make the velocity field divergence-free again

    free(w); free(eta_x); free(eta_y); free(norm_eta_x); free(norm_eta_y);
}

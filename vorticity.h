#ifndef VORTICITY_H
#define VORTICITY_H

void compute_vorticity(
    int N, 
    float *u, 
    float *v, 
    float *w
    );

void compute_eta(
    int N, 
    float *w, 
    float *eta_x, 
    float *eta_y
    );

void normalize_eta(
    int N, 
    float *eta_x, 
    float *eta_y, 
    float *norm_eta_x, 
    float *norm_eta_y
    );

void vorticity_confinement(
    int N, 
    float *u, 
    float *v, 
    float *w, 
    float *eta_x, 
    float *eta_y, 
    float *norm_eta_x, 
    float *norm_eta_y, 
    float dt, 
    float epsilon
    );


#endif // VORTICITY_H

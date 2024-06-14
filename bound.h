#ifndef BOUND_H
#define BOUND_H

typedef struct {
    int *fixed_cells; // Fixed cells
} BoundaryCells;

// Initializes fixed cells
void init_boundary_cells(BoundaryCells *cells, int N);

// Sets a fixed cell
void set_fixed_cell(BoundaryCells *cells, int N, int i, int j);

// Frees fixed cells
void free_boundary_cells(BoundaryCells *cells);

#endif // BOUND_H

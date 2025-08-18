# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

**Build the demo:**
```bash
make
```

**Run the demo:**
```bash
./demo [N] [dt] [diff] [visc] [force] [source]
# Example: ./demo 128 0.1 0 0 5 100
```

**Clean build artifacts:**
```bash
make clean
```

## Code Architecture

This is a 2D fluid simulation implementation based on Jos Stam's "Stable Fluids" algorithm, enhanced with vorticity confinement and fixed boundary cells.

### Core Components

**Main solver (`solver.c`):**
- `vel_step(N, u, v, u_prev, v_prev, visc, dt, epsilon, &cells)` - Velocity field evolution with optional vorticity confinement
- `dens_step(N, dens, dens_prev, u, v, diff, dt, &cells)` - Density field evolution
- Uses `(N+2)×(N+2)` grids with boundary cells at edges
- Coordinate conversion via `IX(i,j)` macro: `((i)+(N+2)*(j))`

**Vorticity system (`vorticity.c/h`):**
- `compute_vorticity()` - Calculates curl of velocity field
- `vorticity_confinement()` - Adds rotational forces to maintain fluid detail
- Controlled by `epsilon` parameter in `vel_step()`

**Boundary handling (`bound.h`):**
- `BoundaryCells` struct manages fixed obstacle cells
- `set_fixed_cell()` - Mark cells as solid boundaries
- Automatically enforced in solver steps

**Demo interface (`demo.c`):**
- OpenGL/GLUT visualization
- Mouse interaction: left click adds velocity, right click adds density
- Keyboard controls: `v` velocity view, `c` clear, `1` toggle vorticity, `2` toggle fixed cells

### Key Data Structures

All fluid quantities stored as 1D arrays sized `(N+2)*(N+2)`:
- `u, v` - Current velocity components
- `u_prev, v_prev` - Previous velocity (for sources/forces)
- `dens, dens_prev` - Current and previous density
- Grid indexing: `IX(i,j)` where `i,j ∈ [0, N+1]`, simulation domain is `[1, N]`

### Integration Notes

The solver operates on a regular grid with periodic time stepping. For UE5 integration, the README.md provides detailed guidance on wrapping the C solver in Unreal components, mapping world coordinates to grid space, and connecting to Niagara particle systems.

Both `code/` directory (original) and root directory contain solver implementations - the root version includes enhancements for vorticity confinement and fixed boundaries.
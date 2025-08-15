### Fluid Solver (Stable Fluids) + UE5 Swimming Simulator

- **Repo contents**: 2D fluid solver in C (`solver.c`, `vorticity.c/h`, `bound.h`) and an OpenGL/GLUT demo (`demo.c`, `Makefile`).
- **Key entry points**:
  - `vel_step(N, u, v, u_prev, v_prev, visc, dt, epsilon, &cells)`
  - `dens_step(N, dens, dens_prev, u, v, diff, dt, &cells)`
  - Optional: vorticity confinement via `epsilon`; fixed cells via `BoundaryCells`.

### Build the C demo (Linux)

```bash
cd /workspace
make
./demo
```

Args (optional): `N dt diff visc force source`. Example: `./demo 128 0.1 0 0 5 100`.

### UE5 integration quickstart (swimming)

1. **Create project**: UE5 C++ Third Person project.
2. **Add plugin**: New plugin `FluidSim` (Runtime). Add this repo as a `ThirdParty` folder or copy `solver.c`, `vorticity.c/h`, `bound.h` into the plugin’s `Source/ThirdParty/StableFluids`.
3. **Build rules**: In the plugin `Build.cs`, include the C sources (Compile as C) or link a static lib; expose headers to the module.
4. **Component wrapper**: Implement `UFluidSolverComponent` that owns arrays: `u, v, u_prev, v_prev, dens, dens_prev` sized `(N+2)^2`. On `BeginPlay` allocate/clear; on `TickComponent` call `vel_step` then `dens_step` each frame.
5. **World↔Grid mapping**: Define a water region in world space (plane or box). Map actor positions `(X,Y)` to grid `(i,j)`; scale forces by `dt`.
6. **Character interaction**: When the player swims (inputs or animation notifies), write impulses into `u_prev/v_prev` near the hands/feet; write dye into `dens_prev` for visual feedback.
7. **Forces on character**: Sample local fluid velocity to compute drag; apply buoyancy based on immersion depth; add lateral flow from `(u,v)` to the physics body.
8. **Visuals**:
   - Niagara: Upload `(u,v)` as a RFloat2 texture each tick; Niagara GPU sim samples it as a vector field for water particles.
   - Material: Drive a water plane material with flow maps from `(u,v)`; use `dens` for foam/ink.
9. **Scenarios**:
   - Pond: Closed boundaries (default). Optionally mark obstacles via `BoundaryCells`.
   - River: Apply constant inflow at left boundary and outflow at right; inject velocity sources each tick.
   - Beach: Combine a simple wave displacement (Gerstner/sine) for surface mesh with this 2D flow for near‑shore currents and foam.
10. **Performance tips**: Start with `N=128`, `dt≈0.016–0.033`, `visc/diff=0`. Run solver on a worker thread; upload textures with `RHIUpdateTexture2D`/`UTexture2D::UpdateResource`.

### Controls (demo)
- Left mouse: add velocity. Right mouse: add density. Keys: `v` toggle velocity view, `c` clear, `q` quit. `1` toggle vorticity, `2` toggle fixed cells.
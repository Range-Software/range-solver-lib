# Range FEA Solver Library

Finite element solver library for Range FEA. Provides physics solvers, linear system solvers, and supporting utilities used by the `fea-solver` CLI and the GUI application.

## Solvers

Each physics solver derives from `RSolverGeneric` and operates on an `RModel` from `range-model-lib`.

| Class | Physics |
|---|---|
| `RSolverHeat` | Steady-state and transient heat conduction |
| `RSolverRadiativeHeat` | Radiative heat transfer (hemicube method) |
| `RSolverFluid` | Incompressible fluid flow (Navier-Stokes) |
| `RSolverFluidHeat` | Coupled fluid flow and heat transfer |
| `RSolverFluidParticle` | Particle transport in fluid |
| `RSolverStress` | Linear structural stress analysis |
| `RSolverAcoustic` | Acoustic pressure wave propagation |
| `RSolverElectrostatics` | Electrostatic field analysis |
| `RSolverMagnetostatics` | Magnetostatic field analysis |
| `RSolverWave` | General wave equation |
| `RSolverMesh` | Mesh deformation / moving mesh |

The top-level `RSolver` class owns the solver map and orchestrates multi-physics execution.

## Linear System Solvers

`RMatrixSolver` wraps iterative solvers for sparse systems assembled by the physics solvers:

- Conjugate Gradient (CG)
- Generalized Minimal Residual (GMRES)

`RMatrixPreconditioner` provides preconditioning, and `RMatrixManager` handles sparse matrix storage.

`REigenValueSolver` provides eigenvalue/eigenvector computation for modal analysis.

## Supporting Utilities

- `RConvection` — convective boundary condition assembly
- `RLocalRotation` — local coordinate frame rotations for anisotropic material properties
- `RScales` — physical unit scaling factors
- `RIterationInfo` / `RIterationInfoValue` — convergence tracking per iteration
- `RHemiCube` / `RHemiCubePixel` / `RHemiCubeSector` — hemicube geometry for view-factor computation in radiative heat transfer
- `RSolverSharedData` — data shared between coupled solvers across time steps

## Dependencies

- `range-base-lib` — math primitives, logging, job management
- `range-model-lib` — FEA model: elements, nodes, materials, boundary conditions

## Version 1.2.0

### Improvements

- **RSolverAcoustic:** the solver is functional again and is no longer flagged
  as *NOT WORKING*.
  - Added a frequency domain (harmonic) analysis driven by **RAcousticSetup**.
    The complex Helmholtz system is solved as an equivalent real block system
    of twice the size, producing one result record per swept frequency.
  - The absorbing boundary is now assembled as a boundary damping matrix
    derived from the boundary condition absorption coefficient instead of being
    patched onto the solution after the solve, and a new **Acoustic impedance**
    boundary condition assembles the same term from a specific impedance.
  - Bulk attenuation is supported through the new **Acoustic damping factor**
    material property, and the speed of sound may now be given directly instead
    of being derived from the modulus of elasticity and the density.
  - Results were extended with the sound pressure level, the acoustic
    intensity, the acoustic phase and the imaginary part of the velocity
    potential.
  - Point entities now contribute their lumped boundary terms; previously they
    were silently skipped because point elements carry no integration points.
  - Density and either a speed of sound or a modulus of elasticity are enough
    to make an element computable; the two stiffness sources are alternatives,
    not both required.
- **RSolver::run():** a harmonic acoustic analysis is no longer driven by the
  time loop. The time solver setup is left untouched, so switching back to a
  transient analysis does not lose it.
- **RSolverGeneric::run():** added a frequency sweep branch for harmonic
  acoustics, mirroring the existing modal analysis branch. `writeResults()`
  writes one record per frequency, independent of the time solver output
  frequency.
- **RSolverGeneric::findComputableElements()** is now virtual so that a solver
  can define its own rule for which material properties an element needs.
- Added **tst_solver_acoustic**, which verifies the solver against the closed
  form plane wave and standing wave solutions of a one dimensional duct.
- Added **doc/acoustic_theory_manual.md**, covering the formulation, the
  acoustic parts of the user interface and two worked tutorials.
- **RSolverStress:** the individual stress components are stored as results, in
  global coordinates for volume elements and in the local element frame for
  surface and line elements.
- **RSolverGeneric::generateVariableVector():** a condition component which is
  switched off no longer prescribes a value. This makes the optional components
  of the *Displacement* boundary condition work, so a support can hold one
  global direction and leave the others free.
- Added **tst_solver_stress** and **tst_eigen_value_solver**, which verify the
  static bar solution, the von Mises invariant against the stored stress
  components, the axial modes of a fixed-free bar against `c/(4*L)`, the eigen
  values of a small system against their closed form, that an entered local
  direction changes which degree of freedom a roller restrains, that a
  roller adds a reaction only along the direction it restrains, that constraints
  from two entities combine on a shared node and carry their prescribed value
  into the node frame, and that contradicting constraints are reported.
- **REigenValueSolver:** the multiple mode method was replaced by **subspace
  iteration** with a Rayleigh-Ritz projection, and the dominant mode method by
  an **inverse power iteration** with a Rayleigh quotient. Both converge towards
  the lowest eigen values and return `lambda` of `K*phi = lambda*M*phi`
  directly. The projected eigenproblem is reduced with a Cholesky decomposition
  of the projected mass matrix and solved with the cyclic Jacobi method. The
  `Arnoldi` and `Rayleigh` methods of **REigenValueSolverConf** were renamed to
  `SubspaceIteration` and `InversePowerIteration` accordingly.
- **RSolverStress:** displacement constraints are now resolved per node instead
  of per boundary condition. **RSolverStress::generateLocalConstraints()**
  collects every constraint acting on a node - each of them a statement
  `d . u = v` about one direction - and reduces the collection to at most three
  mutually perpendicular held directions by Gram-Schmidt, carrying the
  prescribed values through the same operations. The held directions become the
  leading axes of the node frame, held in the new `nodeConstrainedDirections`
  and `nodePrescribedDisplacement` members, and whatever is left completes the
  frame and stays free. Consequences:
  - Constraints from different entities meeting at one node combine. A face may
    be held in a global direction by *Displacement* and rolled on a tilted plane
    by *Roller displacement* at the same time; previously the later condition
    overwrote the frame of the earlier one, and a globally phrased component was
    then read in somebody else's local frame.
  - A prescribed value is applied in the frame it was given in, so a
    *Normal displacement* or a *Roller displacement* with a non-zero value now
    moves the node by that amount along its own direction.
  - Two entities which prescribe different values in the same direction are
    reported as a conflict naming the node, instead of being resolved in favour
    of whichever was read last.
  - **RSolverGeneric::updateLocalRotations()** became virtual; the structural
    solver overrides it, because its frames come from the constraints rather
    than from the geometry of a single boundary condition.

### Bug fixes

- **RSolverStress::assemblyMatrix():** the mass matrix of a modal analysis was
  assembled from the element matrix before the local rotations were applied,
  while the stiffness matrix was assembled after them. A mode shape of a model
  carrying a local frame was therefore computed from a mismatched pair. The
  rotated mass matrix is now used.

- **RSolverStress:** the von Mises stress was reported as `QN + QS`, the sum of
  the normal and the shear invariant, instead of `sqrt(QN^2 + QS^2)`. The
  reported value was up to about 41 % too high. The two-dimensional shear
  invariant of a surface element also kept the sign of the shear stress.
- **RSolverStress:** the modal setup now holds the natural frequency in Hz,
  converted as `sqrt(lambda)/(2*pi)`, instead of the raw eigen value, which is
  what the 3D view and the report have always labelled it as. A mode the
  iteration could not resolve is reported as zero with a warning in the log.
- **RSolverStress:** the axial stress of a line element was reported as
  `E*A*eps`, which is an axial force, instead of `E*(eps - alpha*dT)`. The
  thermal part of the same expression carried a further factor of the cross
  area.
- **RSolverStress:** nodal forces were recovered as `M*a + K*u`, but the nodal
  acceleration was only ever read back from a stored *Acceleration* result which
  no solver writes and no condition supplies. The inert term and the dead
  acceleration plumbing were removed - the reported force is the internal
  elastic force `K*u`.
- **RSolverStress::generateNodeBook():** a *Displacement* boundary condition
  constrained all three degrees of freedom regardless of which of its components
  were switched on.
- **REigenValueSolver::solve():** the inversion and ascending sort that turn the
  raw iteration values into the eigen values of `K*phi = lambda*M*phi` were
  applied only when more than one value was extracted, so a single extracted
  value came back on the reciprocal scale. Every method now returns `lambda`
  directly and `solve()` only sorts.
- **REigenValueSolver:** the extracted eigen values were wrong. On a two degree
  of freedom system with known eigen values the old Arnoldi and QR iteration was
  off by tens of percent and varied from run to run, and on a fixed-free bar the
  reported fundamental was closer to the second mode. The dominant mode method
  started its shift at `1e9 * rand()`, which drove it towards the highest modes
  rather than the lowest. Both were replaced, see above. The axial modes of a
  bar now come out within 0.5 % of `(2n+1)*c/(4*L)`.
- **RSolverStress::prepare():** the stiffness matrix of a line element was
  overwritten at every integration point instead of being accumulated, and was
  never multiplied by the Jacobian determinant and the integration weight. A
  truss came out far too stiff - by a factor of the element count on a uniform
  bar. The thermal expansion force of the same element carried an extra factor
  of the cross area and indexed the strain-displacement vector by the node
  number instead of the degree of freedom.
- **RSolverStress::applyLocalRotations():** the element load vector of a node
  carrying a local frame was rotated with the transformation which maps local to
  global, while the element matrix was rotated with its transpose. The solution
  then satisfied equilibrium with a rotated load instead of the applied one, so
  a roller support carrying a load - self weight, a traction, a thermal load -
  came out wrong. It went unnoticed because the rotation matrix is symmetric
  whenever the local direction is along a global axis, which is the usual case.
- **RSolverGeneric::updateLocalRotations():** the local direction entered with a
  *Normal displacement* or a *Roller displacement* was honoured on point
  entities only - a surface always used its averaged element normals and a line
  its element direction, silently ignoring what the user had entered. When the
  boundary condition asks for its direction to be used, it is now applied to
  surfaces and lines as well. A zero length direction is reported as an error
  instead of producing a degenerate frame.
- **RSolverStress::prepare():** *Force* and *Weight* assigned to a point entity
  were applied in full at every one of its point elements. They are totals over
  the entity, as they are for a line or a surface, and are now spread over its
  points.

- **RSolverAcoustic:** the mass matrix was assembled with a negative sign, so
  the transient system matrix `K - a0*M` was indefinite and could not be solved
  by the conjugate gradient method it was handed to.
- **RSolverAcoustic:** the Newmark velocity and acceleration were recovered but
  never stored, so the time integration state was reset to zero at the start of
  every time step. Both are now published as results.
- **RSolverAcoustic:** the Newmark predictor used the boundary condition values
  of the current step instead of the state at the beginning of the step.
- **RSolverAcoustic:** acoustic pressure is evaluated as `p = rho * dphi/dt`;
  it used to be evaluated as `rho * phi`, which is dimensionally wrong.
- **RSolverAcoustic:** the absorbing boundary normal search indexed the edge
  element list with the global element index, reading out of bounds.
- **RSolverAcoustic:** the particle velocity gradient was scaled by the element
  Jacobian determinant and, for volume elements, only the last integration
  point was used.
- **RSolverAcoustic:** the prescribed velocity boundary condition picked up the
  velocity component of unrelated boundary conditions, such as forced
  convection.
- **RSolverAcoustic:** a zero time-march approximation coefficient made the
  Newmark coefficients divide by zero. The coefficient is now clamped to the
  unconditionally stable average acceleration scheme.
- **RSolverAcoustic:** a transient analysis without an enabled time solver used
  to silently degenerate into a Laplace problem; it is now reported as an error.
- **RSolverAcoustic:** the frequency domain solve widens the GMRES restart to
  fit the system, within a fixed memory budget. The default restart of 10 made
  the indefinite block system stagnate.
- **RScales:** the velocity potential and its time derivatives were not scaled
  with the mesh, and the acoustic particle velocity was scaled by `m*s` instead
  of `m/s`.

---

## Version 1.1.0

### Improvements

- Class **RFileManager** changed to namespace **RFileUtils**
- **RSolverAcoustic, RSolverElectrostatics, RSolverHeat, RSolverMagnetostatics,
  RSolverStress:** element assembly no longer runs inside an `omp critical`
  section. Each thread assembles into its own `RSparseMatrix`/`RRVector`
  (plus `M` for modal stress), and the buffers are merged into the global
  system in a single parallel row-wise pass after the element loop. The
  protected `assemblyMatrix()` methods now take the target matrix/vector as
  parameters instead of writing to the solver members directly.
- **RSolverFluid:** per-thread assembly buffers moved into members
  (`threadAssemblyMatrices`, `threadAssemblyVectors`); the sparse pattern is
  copied into them only when the pattern, mesh, or thread count changes, and
  values are zeroed in place between iterations.
- **RSolverAcoustic, RSolverElectrostatics, RSolverFluidHeat,
  RSolverFluidParticle, RSolverHeat, RSolverMagnetostatics, RSolverStress:**
  `bool` abort flags with `#pragma omp flush` replaced by `std::atomic<bool>`
  with `memory_order_relaxed`, removing explicit memory barriers from the
  element loops.
- **RSolverRadiativeHeat::prepare():** matrix rows are pre-allocated with
  `setNRows()`, so each outer iteration touches only its own row and the
  `omp critical` section around `A.addValue()` is gone.
- **RSolverStress::process():** per-element stress results are written outside
  the `omp critical` section that accumulates node forces.
- **RScales::convert():** the condition-component list is gathered once before
  the parallel region instead of being rebuilt by every thread inside an
  `omp critical` section.
- **RSolverMesh::prepare():** maximum element volume uses an OpenMP
  `reduction(max:)` clause instead of an `omp critical` section.
- **RHemiCube::calculateViewFactors():** the eye-patch loop uses
  `schedule(dynamic)`; emitter patches do far more work than non-emitters, so a
  static split left threads idle.
- **RMatrixSolver:** the sparse CSR cache is now guarded by a mutex and bounded
  to 16 entries (entries are keyed by matrix address, so stale ones would
  otherwise accumulate). `solveCG()`/`solveGMRES()` reuse the cache without
  re-validating it, as `solve()` refreshes it immediately beforehand.
- **RMatrixSolver::solveCG(), solveGMRES():** iterations now also stop when the
  solution diverges, not only when it converges.
- **RIterationInfo:** new `hasDiverged()` method reporting a non-finite error or
  trend.
- **RSolverAcoustic, RSolverElectrostatics, RSolverFluid, RSolverFluidHeat,
  RSolverFluidParticle, RSolverHeat, RSolverMagnetostatics,
  RSolverRadiativeHeat:** `catch (RError error) { ...; throw error; }` replaced
  with `catch (const RError &) { ...; throw; }`, avoiding an exception copy and
  preserving the original exception.
- **RSolverFluid::statistics():** removed temporary `VALIDATION` norm logging.

### Bug fixes

- **RIterationInfo::hasConverged():** a non-finite (NaN/infinite) error or trend
  was reported as *converged*, silently accepting a diverged solve. Non-finite
  values now report divergence via `hasDiverged()`.
- **RConvection::calculateNu():** the Churchill & Chu laminar branch for vertical
  planes and cylinders tested `Ra <= 1.0e-9` instead of `Ra <= 1.0e9`, so
  practically every case took the turbulent correlation.
- **RConvection::calculateNu():** the horizontal-plates case dropped the
  unreachable `0.27 * Ra^(1/4)` branch; that correlation needs plate-orientation
  information which is not available here.
- **REigenValueSolver::solve():** eigenvectors were reordered by pairwise swaps
  (preceded by a spurious `d[0]`/`d[1]` swap), which did not reproduce the sort
  permutation. Rows are now permuted through `indexes` into a copy, guarded by a
  dimension check.
- **REigenValueSolver::solveRayleigh():** the shifted system was solved with `M`
  instead of the shifted matrix `M2`.
- **REigenValueSolver::qlDecomposition():** convergence test `m >= l` relaxed the
  exit condition and could terminate the sweep early; corrected to `m == l`.
- **REigenValueSolver::qrDecomposition():** `RLogger::unindent()` was skipped on
  the converged path, leaving log indentation unbalanced.
- **RHemiCube::_init():** existing sectors were leaked when copying into an
  already-populated hemicube; they are now deleted first. `operator=()` also
  guards against self-assignment.
- **RSolverHeat::prepare():** the `B` matrix was not zeroed between line elements,
  so contributions accumulated across elements.
- **RSolverHeat::prepare():** `htc`/`htt` were shared across the parallel surface
  loop while `getNaturalConvection()` overwrites them per element — a data race.
  Each iteration now takes its own copy of the surface values.
- **RSolverGeneric::writeResults():** the time-step modulo was evaluated without
  checking the output frequency, dividing by zero when it is 0.
- **RSolverRadiativeHeat::process():** element/patch area ratio was computed
  without checking for a zero patch area.
- **RMatrixSolver::solveCG():** `ro1 = std::max(ro1,eps)` flipped the sign of a
  negative `ro1`; the magnitude is now clamped while the sign is preserved.
- **RSolverFluid::clearShapeDerivatives():** the pointer vector was not cleared
  after deleting its contents, leaving dangling pointers. `prepare()` now clears
  the cached derivatives before recomputing them when the mesh changed.

---

## Version 1.0.1

### Improvements

- Added unit tests based on QTest framework

---

## Version 1.0.0

### Improvements

- **RSolverFluid::prepare():** `Ae`/`be` element matrices lifted into `FluidMatrixContainer`
  as thread-local storage, eliminating a heap allocation and deallocation for every
  element on every assembly pass.
- **RSolverFluid::prepare():** four sequential `convertNodeToElementVector` calls
  replaced with a single OpenMP parallel loop, reducing four element traversals
  to one.
- **RSolverFluid::prepare():** `bool` abort flag with `#pragma omp flush` replaced by
  `std::atomic<bool>` with `memory_order_relaxed`, removing explicit memory barriers
  from the check at the top of every element iteration.
- **RSolverFluid::prepare():** `std::pow(x, 2)` in surface-element gravity-magnitude
  calculation replaced with direct multiplication.
- **RSolverGeneric::findInwardElements():** O(`N_surfaces` × `N_elements`) nested scan
  replaced with a pre-built node-to-volume-element index constructed in one
  O(`N_elements`) pass, reducing the per-surface-element search from O(`N_elements`)
  to O(`local_node_degree`).
- **RSolverGeneric::processMonitoringPoints():** missing `break` added after the
  containing element is found, stopping the O(`N_elements`) `isInside()` scan as
  soon as the result is known.
- **RHemiCubeSector:** `limitBox` (bounding box of the limit polygon) is now cached at
  construction time instead of being recomputed on every `testVisibility()` call.
- **RHemiCube::calculateViewFactors():** element triangulations are pre-computed once
  before the parallel eye-patch loop instead of being recomputed for every
  (eye-patch, element) pair.
- **RHemiCubeTriangleComp:** sort comparator now uses squared distances, eliminating
  all `sqrt` calls from the O(N log N) triangle sort.
- **RSolverRadiativeHeat::prepare():** per-thread `b[i]` is now accumulated in a local
  variable and the `omp critical` section covers only `A.addValue`, removing the
  serialisation of the entire inner j-loop.
- **RSolverHeat::prepare():** `getTimeSolver().getEnabled()` hoisted out of the
  innermost integration-point loop into a single `const bool`.
- **RSolverStress::prepare(), process():** the compound condition
  `getEnabled() || problemType==MODAL` hoisted out of all inner loops.
- **RConvection::calculateGr():** `pow(x, 2.0)` and `pow(x, 3.0)` replaced with direct
  multiplication.
- **RSolverGeneric:** node-book generation now batches disabled positions and
  rebuilds the compacted book in one linear pass, avoiding repeated full-book
  consolidation while applying boundary-condition constraints.
- **RSolverFluid:** node-book generation uses the batched `RSolverGeneric` path, and
  average stream-velocity calculation now uses direct multiplication plus an
  OpenMP reduction instead of `pow()` and atomic accumulation.
- **RSolverFluidHeat and RSolverFluidParticle:** matrix assembly now uses per-thread
  sparse matrices/vectors and merges after the parallel element loop, removing
  the OpenMP critical section around element assembly.
- **RSolverFluidHeat and RSolverFluidParticle:** shape derivatives, node books, and
  stream velocity are reused across iterations unless the mesh or first-iteration
  state requires recomputation.
- **RSolverFluidHeat and RSolverFluidParticle:** per-element velocity-divergence
  temporaries are reused from matrix containers instead of being allocated in
  every element integration path.
- **RSolverFluid::computeElementConstantDerivative():** intermediate element-level
  matrices/vectors are no longer materialized and cleared for every element.
  The final `Ae`/`be` contributions are accumulated directly, reducing memory
  traffic in the constant-derivative fluid assembly path.
- **RSolverFluid::prepare():** the global sparse matrix pattern is now built once
  when the mesh/node book changes and numeric values are zeroed in-place between
  iterations, avoiding repeated sparse insertions for stable systems.
- **RSolverFluid::prepare():** each assembled element now caches its active local
  DOFs and global sparse row/column positions; threaded assembly iterates only
  active entries and writes directly to cached sparse positions.
- **RMatrixSolver:** CG and GMRES cache sparse row indexes and values once per solve,
  avoiding repeated sparse-matrix lookups in matrix-vector products.
- **RMatrixSolver:** CG updates/residual norms and GMRES Arnoldi orthogonalization
  work are now parallelized more broadly, and hot-loop `pow(x, 2)` calls were
  replaced with direct multiplication.
- **RMatrixSolver:** sparse matrix data is now cached in a flat CSR-style layout and
  reused across solver instances for the same `RSparseMatrix` when the row pattern
  is unchanged; CG/GMRES SpMV and matrix norm evaluation use this cache.
- **RMatrixSolver::solveGMRES():** Arnoldi orthogonalization batches all current
  Krylov dot products and applies the correction in a single vector pass,
  reducing synchronization inside the inner iteration.
- **RMatrixPreconditioner:** Jacobi stores inverse diagonal values and applies them
  in parallel; block Jacobi now precomputes inverse blocks once and applies them
  as dense block-vector products instead of solving each block every iteration.
- Public solver headers and implementations now use empty parameter lists `()`
  instead of `(void)` for no-argument functions.

### Bug fixes

- **RConvection::getFluidTemp():** was returning surface temperature (Ts) instead of
  fluid temperature (Tf).
- **RConvection::calculateNu():** laminar internal forced convection Nusselt number
  was dividing by `(Re*Pr)^(2/3)` instead of multiplying by `(Re*Pr)^(1/3)`.
- **RConvection::calculateNu():** unreachable branch in the horizontal-plates case
  corrected; Ra ranges now cover all cases without overlap or gap.
- **RHemiCubeSector::rayTraceTriangle():** depth-skip guard was testing pixel row `i`
  instead of the current pixel (`pixelId`), causing incorrect occlusion decisions.
- **REigenValueSolver::qlDecomposition():** inner Givens rotation loop incremented `i`
  instead of decrementing it, producing an infinite loop or out-of-bounds access.
  Loop variable changed to `signed int` to handle the `l=0` boundary safely.
- **REigenValueSolver::solveRayleigh():** initial shift `mu` was computed from an
  uninitialised vector element; replaced with a bounded random value.
- **R_EIS_PYTHAG macro:** arguments were not parenthesised, causing wrong results
  when expressions with lower precedence than `*` were passed.
- **RSolverRadiativeHeat::prepare():** RHS vector was assembled into `b[j]` (receiving
  patch) instead of `b[i]` (emitting patch), producing a wrong system. The write to
  `b[i]` for the ambient term was also outside the critical section, causing a data
  race under OpenMP.
- **RSolverGeneric::run():** `static local updateScalesDone` persisted for the process
  lifetime, preventing `updateScales()` from being called on subsequent solves.
  Replaced with the existing `firstRun` member.
- **RSolverWave::prepare():** `static local firstTime` was never reset, so initial
  conditions were skipped on every run after the first. Replaced with `firstRun`.
- **RSolverFluid, RSolverFluidHeat, RSolverFluidParticle:** static locals `counter` and
  `oldResidual` in `statistics()` were never reset between solver runs, corrupting
  convergence history output. Promoted to member variables.
- **RSolverFluidHeat:** element heat accumulation now sizes node accumulation arrays
  by node count instead of element count.
- **RSolverFluidParticle:** particle-rate recovery now uses element count for
  element-applied values instead of node count.

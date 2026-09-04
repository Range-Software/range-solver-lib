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

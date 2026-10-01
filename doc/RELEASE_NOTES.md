## Version 1.2.0

### Improvements

#### Fluid solver convergence and iteration

- **RSolverFluid** now reports convergence; `hasConverged()` used to always
  answer no, so a task group with a flow task ran its full iteration count.
  Both of these must hold:
  - the relative Newton increment `||dv|| / ||v||` and `||dp|| / ||p||`
    (logged as `Convergence-V` / `Convergence-P`) is below the task group
    convergence value
  - the residual has fallen to a tenth of its value at the first pass of the
    same solve (each time step is a solve). This catches nearly singular
    systems - no pressure reference, or a mesh too coarse for the Reynolds
    number - which take tiny steps without converging. The tenth is fixed, so
    a loose convergence value cannot disable it; a slowly converging model may
    run its full iteration count, as before. Ratio and target are logged
- The increments replace the signed difference of field norms, which was
  neither positive nor dimensionless. Both norms are taken inside the
  downscaling bracket, so no scale factor is needed and the separate pass over
  the nodes is gone. The convergence graph now falls monotonically
- At least two passes are taken before convergence is reported, so a field
  that has not moved yet (no inflow, first assembly) is not taken as converged.
  The log prints the group target and marks the iteration that reached it
- **RSolverGeneric::run()** takes the task group convergence value as a third
  argument, passed down by **RSolver::runProblemTask()**
- **RSolverFluid** damps its step as `x += omega * dx`, since its matrix is not
  the exact derivative of the residual. `omega` starts at `1` each solve, is
  halved (floor `0.1`) after a pass that raises the residual by more than
  2 % and grows by a quarter after one that lowers it; it is logged as
  `Relaxation`. A run that never overshoots behaves as before. The 2 %
  tolerance avoids stalling on the normal pass-to-pass wander (without it, and
  a `0.01` floor, a run settled at residual ratio `0.48` instead of `0.17`).
  On the reference transient model the residual of the first step, which used
  to climb `2.79`, `3.05`, `3.07`, now retreats to `omega = 0.5` on the second
  pass and descends to `2.48` by the tenth (`3.07` before)
- The residual is computed once per pass in `updateResidualAndRelaxation()` and
  reused by the convergence test and statistics
- `Tsupg`, `Tlsic` and the element length along the flow are evaluated at the
  field being solved rather than at the previous time step velocity, so in
  transient runs the stabilisation no longer lags behind the flow

#### Fluid Jacobian verification

- **RSolverFluid::verifyJacobian()**, enabled by `--verify-jacobian` on the
  solver command line, compares the assembled matrix with a central finite
  difference of the residual (the exact Newton matrix is `-db/dx`, since the
  solver applies `x += dx` from `A*dx = b`). It runs once and reassembles the
  system twice per unknown, so it is meant for meshes of a few elements
- It reports, per block (velocity/velocity, velocity/pressure,
  pressure/velocity, pressure/pressure), the min, median and max ratio of
  assembled to measured entries, grouped so that distinct factors and their
  entry counts stand out, plus the largest disagreement with its node and
  component. A ratio of `1` throughout means the block is correct
- During the sweep the stabilisation parameters, element length, pass counter
  and mesh-changed flag are frozen (the matrix carries no derivative of them,
  and reassembly must not rebuild or resize the system), while the nodal
  acceleration is refreshed for each perturbed assembly (`prepare()` does not
  recompute it, which would drop the whole time-derivative term)
- The difference step is about the cube root of machine epsilon. At `1.0e-7`
  cancellation made the pressure/pressure block of a model at rest read `1.025`;
  it now reads `1`

#### Fluid Jacobian defects - diagnosed, not applied

The check shows four fluid matrix terms that are not the derivative of the
residual. Each correction was written, verified with the check and measured;
**none is applied**, because each makes the iteration worse in practice:

- **Viscous coefficient:** the matrix carries `ro * u`, the residual `u`. On a
  steady model at rest the velocity/velocity block reads exactly `1e+06` (the
  density in solver units) - for water the matrix is a million times too stiff
  in the viscous directions, which is why steady-state models crawl
- **Mass terms:** the matrix uses a consistent mass `ro*iNiN` and consistent
  SUPG mass, the residual a lumped one (`mvScale*ax[m]`, `ctvScale*ax[m]`)
- **PSPG acceleration term:** the same mismatch in `bteScale`
- **Theta weighting:** the matrix scales spatial terms by `alpha*dt`, the
  residual by `dt`, so under the default central difference march the
  velocity/velocity block reads `0.5`

With all four corrected every block reads `1` (bar convective linearisation
scatter), yet:

| configuration | steady-state | reference transient step, 10 passes |
|---|---|---|
| **as shipped** | monotone, no rise in 498 passes | `3.03` |
| viscous and mass corrected | **oscillates - 119 rises in 258 passes**, residual up 75 % on the second pass | `0.557` |
| plus the theta correction | - | `1.96` |
| plus the PSPG mass correction | - | `155`, climbing |

The corrected Jacobian is much softer and the first Newton step overshoots by
three quarters. Slower relaxation growth (`1.05`, `1.02`) made it worse, and a
backtracking line search cost more passes than it saved. The system is
convection dominated with an indefinite true Jacobian, so the wrong terms act
as stabilisation: **a more exact Jacobian is not automatically a better
iteration matrix here.** Applying the corrections needs a residual line search
cheap enough to try several step lengths per pass, which the current structure
cannot do since the residual is only available from a full `prepare()`.

Separately, the general (non-tetrahedral) path writes its mass as
`me[m][m] = ro*N[m]*N[n]` inside the `n` loop, so the last `n` wins: hexahedral
transient models assemble a meaningless mass matrix. Tetrahedral models are
unaffected.

#### Magnetostatic solver

- **RSolverMagnetostatics** evaluates the field by the **Biot-Savart law**
  instead of assembling `laplace(B) = -u0 * curl(J)`, which had no boundary
  condition, so every model was singular and the field level was set by the
  iteration rather than the physics. There is now no system to solve:
  `assemblyMatrix()` and the matrix solver setup are gone; `RSolverGeneric::u0`
  is used as before
- The per-element constant current density is integrated in closed form:
  tetrahedra reduce to their faces by the gradient theorem (edge logarithms
  plus face solid angle), surface triangles are sheet currents
  `J * thickness`, quadrilaterals are split into two triangles, and two-node
  segments are straight wires carrying `J * cross area`. Beyond twice the
  longest element edge a degree-two rule is used (3 points on a triangle, 4 on
  a tetrahedron), beyond six times the midpoint rule; the error summed over a
  conductor stays below `1e-4` of the peak field at every node
- Nodes on a volume conductor surface get the (continuous) surface field. On a
  line or surface conductor, where the field is singular, a regularised value
  is returned: a segment contributes nothing on its own line, and a sheet's
  edge logarithm is omitted at nodes on that edge
- **Surface and line entities now carry current** (previously volumes only),
  and **every node receives a field**, including nodes of surface, line and
  point entities (previously unknowns with a zero diagonal) and of a mesh
  around the conductor carrying no current
- Current density is used per element as stored instead of being averaged to
  nodes and interpolated back, which smeared it across conductor edges.
  Elements below `1e-10` of the largest current density (e.g. air beside
  copper) are skipped; with no current the log says so and the field is zero
- Cost is nodes times current-carrying elements; nodes are processed in
  parallel blocks of 64 so each source element is read once per block. On
  `120 000` tetrahedra and `24 000` nodes the evaluation dropped from about 6.4
  to 2.2 s (midpoint rules only) and the full run takes about 5 s on 14 threads
- Only the modelled current contributes: the leads closing an open current path
  between electrodes of an electrostatic model are not included
- Unit tests in `tst_solver_magnetostatics` check the closed forms against
  high-order Gauss-Legendre quadrature, field continuity at tetrahedron
  vertices, edges and faces, the jump across a current sheet, the quadrature
  rules against the closed form over a bar, a bar near and far against direct
  integration and the finite wire formula, a wire, a strip, and an
  electrostatics task driving a magnetostatics one

#### Conjugate heat transfer at fluid walls

- *Forced convection* on a wall between a meshed solid and a meshed fluid now
  couples **RSolverHeat** and **RSolverFluidHeat** directly instead of feeding
  the flat plate correlation an average of the fluid element behind the wall,
  most of whose nodes lie on the wall: the velocity came out as a fraction of
  the off-wall node (zero where it sat on another wall, leaving the wall
  unconvected) and the fluid temperature was pulled towards the wall's
- **RSolverFluidHeat** holds the wall nodes at the solid temperature from the
  heat solver (insulated before the first heat solve) and publishes per wall
  element a heat transfer coefficient `k * G` - the first fluid element
  conductance, `G` being the sum of the off-wall node shape function
  derivatives along the wall normal (the reciprocal element height for a
  tetrahedron) - and a reference temperature reproducing the heat entering the
  fluid at the wall nodes, taken from the fluid system residual. Using the first
  element gradient instead missed the unit test interface temperature by 1.5 K.
  Before the first heat solve the reference temperature is the
  gradient-weighted temperature of the off-wall nodes
- **RSolverHeat** applies the pair as a *Simple convection* condition. The
  correlation, with the configured *Fluid temperature* and *Velocity*, is used
  only on elements no fluid heat result covers (no meshed fluid, or the first
  coupled pass); the log reports per entity which is in use and on how many
  elements
- The wall temperature passed to the fluid is **Aitken** relaxed, bounded to
  `[-100, 100]` and restarted each task run. Plain alternation contracts by
  `(h - S) / (Ks + h)` per pass (`h` first element, `S` whole fluid, `Ks` solid
  conductance), close to one for a poor conductor against a resolved fluid:
  the unit test needs about 160 passes without relaxation and 5 with it
- While coupled, both solvers report `||dT|| / ||T||` as convergence; otherwise
  they stay unconditionally converged, so a group with only an uncoupled heat
  task ends after one iteration. **RSolverFluidHeat** writes the relative change
  to its convergence file instead of the difference of field norms
- Coupling data is shared in SI units under
  `RSolverFluidHeat::wallHeatTransferCoefficientKey`,
  `RSolverFluidHeat::wallFluidTemperatureKey` and
  `RSolverHeat::solidNodeTemperatureKey`. Wall pairing moved to
  **RSolverFluidHeat::findWallElements()** and covers only walls with *Forced
  convection*; `RSolverFluidHeat::fluidNodeTemperatureKey`,
  `fluidNodeVelocityKey`, `RSolverHeat::findFluidTemperature()`,
  `findFluidVelocity()` and `findFluidElements()` are gone.
  **RSolverHeat::getForcedConvection()** is now per element
- **RSolverHeat** solves **solids only**: `findComputableElements()` drops every
  volume whose material **RMaterial::isFluid()** reports as fluid, whatever it
  carries (mercury, carrying an emissivity, used to be conducted as a solid,
  overwriting the fluid heat result). Point, line and surface elements made
  computable only by a condition are dropped where all their nodes touch the
  fluid and not all touch the solid (inlets, outlets); a condition on them the
  fluid heat solver does not read is reported
- **RSolverHeat::checkConvectionInput()** new: a correlated convection value
  that leaves the correlation nothing to work with (zero dynamic viscosity,
  thermal conductivity, hydraulic diameter, density, heat capacity or mean
  velocity) stops the solver naming the component and entity, instead of
  yielding `h = 0` through the **RConvection** division guards
- **RSolverGeneric::generateHeatVector()** and **findElementGroupMeasure()** new,
  spreading the *Heat* boundary condition total `[W]` over the volume, area,
  length or point count of the entity's computable elements (see bug fixes)
- Unit tests in `tst_solver_heat_coupling` check that a heat task leaves a fluid
  carrying every heat property untouched, the parabola of a uniformly heated
  fluid at rest, and the interface temperature of a solid slab against a fluid
  slab - at rest and flowing towards the wall at Peclet number 7 - against the
  closed form within ten coupled passes

#### Acoustic solver

- **RSolverAcoustic** is functional again and no longer flagged *NOT WORKING*:
  - new frequency domain (harmonic) analysis driven by **RAcousticSetup**,
    solving the complex Helmholtz system as a real block system of twice the
    size, one result record per swept frequency
  - the absorbing boundary is assembled as a boundary damping matrix from the
    absorption coefficient instead of being patched onto the solution; the new
    **Acoustic impedance** condition assembles the same term from a specific
    impedance
  - bulk attenuation via the new **Acoustic damping factor** material property;
    the speed of sound may be given directly instead of being derived from
    modulus of elasticity and density - density plus either one makes an
    element computable
  - new results: sound pressure level, acoustic intensity, acoustic phase and
    imaginary velocity potential
  - point entities contribute their lumped boundary terms (previously skipped,
    having no integration points)
- **RSolver::run()** no longer drives a harmonic acoustic analysis through the
  time loop and leaves the time solver setup untouched
- **RSolverGeneric::run()** has a frequency sweep branch for harmonic acoustics,
  like the modal branch; `writeResults()` writes one record per frequency
  regardless of the time solver output frequency
- **RSolverGeneric::findComputableElements()** is now virtual, so a solver can
  define which material properties an element needs
- Added **tst_solver_acoustic** (closed form plane and standing waves in a 1D
  duct) and **doc/acoustic_theory_manual.md** (formulation, user interface and
  two tutorials)

#### Stress and modal analysis

- **REigenValueSolver** replaces the multiple mode method by **subspace
  iteration** with Rayleigh-Ritz projection (projected problem reduced by
  Cholesky of the projected mass and solved by cyclic Jacobi), and the dominant
  mode method by **inverse power iteration** with a Rayleigh quotient. Both
  converge to the lowest eigenvalues and return `lambda` of
  `K*phi = lambda*M*phi` directly. **REigenValueSolverConf** `Arnoldi` and
  `Rayleigh` renamed to `SubspaceIteration` and `InversePowerIteration`
- **RSolverStress** resolves displacement constraints per node:
  **generateLocalConstraints()** collects every constraint `d . u = v` on a node
  and reduces them by Gram-Schmidt to at most three perpendicular held
  directions with their prescribed values (new `nodeConstrainedDirections` and
  `nodePrescribedDisplacement`); the remaining axes complete the frame and stay
  free. As a result:
  - constraints from different entities on one node combine (e.g. a global
    *Displacement* plus a tilted *Roller displacement*) instead of the later
    overwriting the earlier frame and global components being read in a
    foreign local frame
  - prescribed values apply in their own frame, so a non-zero *Normal* or
    *Roller displacement* moves the node by that amount along its direction
  - different values prescribed in the same direction are reported as a
    conflict naming the node, instead of the last one read winning
  - **RSolverGeneric::updateLocalRotations()** is virtual and overridden by the
    stress solver, whose frames come from the constraints
- **RSolverStress** stores individual stress components as results - global
  coordinates for volumes, local element frame for surfaces and lines
- **RSolverGeneric::generateVariableVector()**: a switched-off condition
  component no longer prescribes a value, so *Displacement* can hold one global
  direction and leave the others free
- Added **tst_solver_stress** and **tst_eigen_value_solver**, verifying the
  static bar, von Mises against the stored components, fixed-free bar axial
  modes against `c/(4*L)`, eigenvalues of a small system against the closed
  form, that an entered local direction selects the degree of freedom a roller
  restrains and the reaction is only along it, that constraints from two
  entities combine on a shared node with their prescribed values, and that
  contradicting constraints are reported

### Bug fixes

#### Fluid heat

- **RSolverFluidHeat** assembled conduction (`-k*grad(N).grad(N)`) with the
  opposite sign to advection (`+rho*c*v.grad(T)`), i.e. the flow ran
  backwards: heat was carried upstream and a *Heat* source cooled the fluid.
  A fluid at rest without a source was unaffected

#### Electrostatics

- **RSolverElectrostatics** recovered the field gradient in `process()` weighted
  by the element Jacobian determinant (and, on surfaces, the thickness), so
  electric field, current density, electric energy and Joule heat scaled with
  element size and did not converge under refinement. The nodal potential and
  the resistivity `|E|/|J|` were unaffected. Thickness and cross area remain in
  the stiffness, and the `getThickness() > 0` / `getCrossArea() > 0` guards
  still skip non-conducting entities
- **Charge density** is assembled with a positive sign on every element type
  (from `div(e0*er*grad(V)) = -rho`); line, surface and volume loops subtracted
  it, so a positive charge lowered the potential and point and volume charges
  acted oppositely. Models driven only by prescribed potentials were unaffected
- Joule heat is stored as the dissipation density `sigma*|E|^2` in `W/m^3`, as
  **RSolverHeat** and **RSolverFluidHeat** expect; it carried an extra
  characteristic element length, making the power of a resistive heater mesh
  dependent. The length computation is removed, and the (unused) **RScales**
  dimension changed from `kg*m^2/s^3` to `kg/(m*s^3)`

#### Heat

- **RSolverHeat** and **RSolverFluidHeat**: the *Heat* boundary condition,
  labelled `[W]`, acted as `W/m^3`, `W/m^2` or `W/m` on volumes, surfaces and
  lines, and was scaled as a total power, so the delivered power depended on
  model size (`1000 W` on a `0.01 m^2` face delivered `40 W`). The total is now
  spread over the entity and the energy balance closes
- **RSolverHeat::getForcedConvection()** demanded a *Fluid temperature*
  component the condition was created without, so any model using *Forced
  convection* stopped with an error. The component now exists and is the
  fall-back where no fluid heat result covers the wall

#### Stress and modal analysis

- **RSolverStress**:
  - von Mises was `QN + QS` instead of `sqrt(QN^2 + QS^2)` (up to about 41 %
    too high); the 2D shear invariant of surface elements kept the shear sign
  - modal frequencies are stored in Hz (`sqrt(lambda)/(2*pi)`) as the 3D view
    and report label them, instead of the raw eigenvalue; an unresolved mode is
    reported as zero with a warning
  - line element axial stress was `E*A*eps` (a force) instead of
    `E*(eps - alpha*dT)`, with an extra cross area factor on the thermal part
  - nodal forces were `M*a + K*u` with an acceleration nothing ever wrote; the
    dead term and plumbing are removed and the force is `K*u`
  - `generateNodeBook()`: *Displacement* constrained all three degrees of
    freedom regardless of which components were on
  - `prepare()`: line element stiffness was overwritten at each integration
    point and never multiplied by the Jacobian determinant and weight, making
    trusses too stiff (by the element count on a uniform bar); its thermal
    force had an extra cross area factor and indexed the strain-displacement
    vector by node instead of degree of freedom
  - `prepare()`: *Force* and *Weight* on a point entity were applied in full at
    every point element; they are now spread over its points
  - `assemblyMatrix()`: the modal mass matrix was assembled before local
    rotations and the stiffness after; the rotated mass is now used
  - `applyLocalRotations()`: the load vector at a local-frame node was rotated
    with the transpose of the matrix rotation, so loaded rollers (self weight,
    traction, thermal load) were wrong whenever the local direction was not
    along a global axis
- **RSolverGeneric::updateLocalRotations()**: an entered *Normal* or *Roller
  displacement* direction was honoured only on points; surfaces and lines used
  their geometry. It now applies to them when the condition asks for it, and a
  zero length direction is an error
- **REigenValueSolver**: eigenvalues were wrong - the old Arnoldi/QR iteration
  was off by tens of percent and non-deterministic on a two-degree-of-freedom
  system, and the dominant mode method's shift started at `1e9 * rand()`,
  driving it to the highest modes. With the new methods bar axial modes are
  within 0.5 % of `(2n+1)*c/(4*L)`. `solve()` inverted and sorted the values
  only when more than one was extracted (a single value came back reciprocal);
  methods now return `lambda` and `solve()` only sorts

#### Acoustics

- **RSolverAcoustic**:
  - the mass matrix had a negative sign, making `K - a0*M` indefinite and
    unsolvable by conjugate gradients
  - Newmark velocity and acceleration were never stored, resetting the time
    integration state every step; both are now results
  - the Newmark predictor used current step boundary values instead of the
    start-of-step state
  - acoustic pressure is `p = rho * dphi/dt`, not `rho * phi`
  - the absorbing boundary normal search indexed edge elements by global element
    index, reading out of bounds
  - the particle velocity gradient was scaled by the Jacobian determinant and,
    for volumes, used only the last integration point
  - prescribed velocity picked up velocity components of unrelated conditions
    such as forced convection
  - a zero time-march coefficient divided by zero; it is clamped to the
    average acceleration scheme
  - a transient analysis without an enabled time solver silently became a
    Laplace problem; it is now an error
  - the frequency domain GMRES restart (default 10, which stagnated) is widened
    to fit the system within a fixed memory budget
- **RScales**: velocity potential and its time derivatives were not scaled with
  the mesh, and acoustic particle velocity was scaled by `m*s` instead of `m/s`

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

#ifndef RSOLVERACOUSTIC_H
#define RSOLVERACOUSTIC_H

#include "rsolvergeneric.h"

//! Acoustic solver.
//!
//! Solves the (damped) wave equation for the velocity potential phi
//!
//!     c^2 * laplace(phi) - d2phi/dt2 - beta * dphi/dt = 0
//!
//! discretised as  M * phi'' + C * phi' + K * phi = f  with
//!
//!     M = int( N_m * N_n )                 mass
//!     K = c^2 * int( grad(N_m).grad(N_n) ) stiffness
//!     C = beta * M + boundary damping      damping
//!     f = c^2 * int( v_n * N_m )           prescribed normal velocity source
//!
//! Sign conventions: particle velocity u = -grad(phi) and acoustic pressure
//! p = rho * dphi/dt. A prescribed velocity boundary condition is positive when
//! it drives the acoustic domain (a piston pushing into the fluid).
//!
//! In transient mode the system is integrated with the Newmark scheme. In
//! harmonic mode phi = Re{ Phi * exp(i*omega*t) } and the complex system
//! (K - omega^2 * M + i*omega*C) * Phi = F is solved as an equivalent real
//! block system of twice the size.
class RSolverAcoustic : public RSolverGeneric
{

    protected:

        //! Element speed of sound.
        RRVector elementSoundSpeed;
        //! Element density.
        RRVector elementDensity;
        //! Element bulk damping factor.
        RRVector elementDampingFactor;
        //! Element boundary damping coefficient (absorbing / impedance boundary).
        RRVector elementBoundaryDamping;
        //! Element prescribed normal velocity.
        RRVector elementVelocityNormal;

        //! Node velocity potential (real part in harmonic mode).
        RRVector nodeVelocityPotential;
        //! Node velocity potential at the beginning of the current time step.
        RRVector nodeVelocityPotentialOld;
        //! Node velocity potential imaginary part (harmonic mode only).
        RRVector nodeVelocityPotentialImag;
        //! Node velocity potential - first time derivative.
        RRVector nodeVelocityPotentialVelocity;
        //! Node velocity potential - second time derivative.
        RRVector nodeVelocityPotentialAcceleration;

        //! Node acoustic pressure (amplitude in harmonic mode).
        RRVector nodeAcousticPressure;
        //! Node acoustic pressure phase in degrees (harmonic mode only).
        RRVector nodeAcousticPressurePhase;
        //! Node sound pressure level.
        RRVector nodeSoundPressureLevel;
        //! Element acoustic particle velocity (amplitude in harmonic mode).
        RSolverCartesianVector<RRVector> elementAcousticParticleVelocity;
        //! Element acoustic particle velocity imaginary part (harmonic mode only).
        RSolverCartesianVector<RRVector> elementAcousticParticleVelocityImag;
        //! Element acoustic intensity (time averaged in harmonic mode).
        RSolverCartesianVector<RRVector> elementAcousticIntensity;

        //! Indicator whether a frequency-domain analysis is being solved.
        bool harmonic;
        //! Angular frequency of the harmonic analysis [rad/s].
        double angularFrequency;

    public:

        //! Constructor.
        explicit RSolverAcoustic(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverAcoustic();

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        //! Find average speed of sound over all computable elements.
        double findAverageSoundSpeed() const;

        //! Find computable elements.
        //! Acoustics needs a density and either a speed of sound or a modulus
        //! of elasticity - the two are mutually exclusive alternatives, so the
        //! generic "all properties are required" rule cannot be used.
        void findComputableElements(RProblemType problemType) override;

        //! Initialize solver.
        void initialize() override;

        //! Update scales.
        void updateScales() override;

        //! Recover previously computed results.
        void recover() override;

        //! Prepare solver.
        void prepare() override;

        //! Run matrix solver.
        void solve() override;

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Generate element material and boundary condition vectors.
        void generateElementVectors();

        //! Generate element vector from one specific boundary condition type.
        //! Unlike RSolverGeneric::generateVariableVector this does not pick up
        //! same-named components of boundary conditions belonging to other
        //! physics (e.g. the velocity of a forced convection boundary).
        void generateBoundaryConditionVector(RBoundaryConditionType boundaryConditionType,
                                             RVariableType variableType,
                                             RRVector &values,
                                             RBVector &setValues) const;

        //! Find element speed of sound and density, filling elements without
        //! material properties with the computable domain average.
        void generateMaterialVectors();

        //! Find boundary damping coefficient for each element from the
        //! absorbing boundary / acoustic impedance boundary conditions.
        void generateBoundaryDamping();

        //! Assembly element matrices into the global matrix system.
        void assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Assembly element matrices for a transient (Newmark) time step.
        void assemblyMatrixTransient(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Assembly element matrices for a harmonic (frequency domain) solve.
        void assemblyMatrixHarmonic(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Run matrix solver for a transient time step.
        void solveTransient();

        //! Run matrix solver for a harmonic frequency.
        void solveHarmonic();

        //! Process acoustic pressure, phase and sound pressure level.
        void processAcousticPressure();

        //! Process acoustic particle velocity.
        void processAcousticParticleVelocity();

        //! Process acoustic intensity.
        void processAcousticIntensity();

        //! Return Newmark gamma coefficient.
        static double findNewmarkGamma();

        //! Return Newmark beta coefficient for given time-march approximation.
        double findNewmarkBeta() const;

        //! Return current time step size (guaranteed to be positive).
        double findTimeStepSize() const;

};

#endif // RSOLVERACOUSTIC_H

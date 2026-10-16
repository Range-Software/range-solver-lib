#ifndef RSOLVERSTRESS_H
#define RSOLVERSTRESS_H

#include <rbl_rvector.h>

#include "rsolvergeneric.h"

class RSolverStress : public RSolverGeneric
{

    protected:

        //! Per-thread assembly buffers.
        struct AssemblyBuffer
        {
            //! System matrix.
            RSparseMatrix A;
            //! Mass matrix (modal analysis and initial acceleration).
            RSparseMatrix M;
            //! Acceleration part of the Newmark predictor (initial acceleration).
            RSparseMatrix D;
            //! Right hand side.
            RRVector b;
            //! Initial equilibrium residual f - K*u - C*v (initial acceleration).
            RRVector r;
            //! Node load in global frame, 3 components per node.
            RRVector f;
        };

        //! Element modulus of elasticity vector.
        RRVector elementElasticityModulus;
        //! Element Poisson ratio.
        RRVector elementPoissonRatio;
        //! Element density vector.
        RRVector elementDensity;
        //! Element thermal expansion coefficient vector.
        RRVector elementThermalExpansion;
        //! Element Rayleigh mass damping coefficient vector.
        RRVector elementRayleighMass;
        //! Element Rayleigh stiffness damping coefficient vector.
        RRVector elementRayleighStiffness;
        //! Element environment temperature vector.
        RRVector elementEnvironmentTemperature;
        //! Number of constrained directions of each node, 0 to 3. They are the
        //! first directions of the node local frame held in localRotations.
        RUVector nodeConstrainedDirections;
        //! Prescribed displacement of each node, expressed in the frame of that
        //! node - the local frame where one is active, global otherwise. Only
        //! the first nodeConstrainedDirections components carry a meaning.
        RRMatrix nodePrescribedDisplacement;
        //! Node displacement vector.
        RSolverCartesianVector<RRVector> nodeDisplacement;
        //! Node initial displacement vector.
        RSolverCartesianVector<RRVector> nodeInitialDisplacement;
        //! Node velocity vector.
        RSolverCartesianVector<RRVector> nodeVelocity;
        //! Node acceleration vector.
        RSolverCartesianVector<RRVector> nodeAcceleration;
        //! Node displacement at the beginning of the time step.
        RSolverCartesianVector<RRVector> nodeDisplacementOld;
        //! Node velocity at the beginning of the time step.
        RSolverCartesianVector<RRVector> nodeVelocityOld;
        //! Node acceleration at the beginning of the time step.
        RSolverCartesianVector<RRVector> nodeAccelerationOld;
        //! Node load of the current time step in global frame.
        RSolverCartesianVector<RRVector> nodeLoad;
        //! Node load of the previous time step in global frame.
        RSolverCartesianVector<RRVector> nodeLoadOld;
        //! True if nodeLoadOld holds the load of the previous time step.
        bool nodeLoadOldValid;
        //! Time step whose beginning is held in the *Old vectors.
        uint oldStateTimeStep;
        //! Node force vector.
        RSolverCartesianVector<RRVector> nodeForce;
        //! Node pressure.
        RRVector nodePressure;
        //! Element stress components.
        //! Volume elements store them in global coordinates, surface and line
        //! elements in their own local element frame.
        //! Order: xx, yy, zz, yz, xz, xy.
        RRVector elementStress[6];
        //! Element normal stress.
        RRVector elementNormalStress;
        //! Element shear stress.
        RRVector elementShearStress;
        //! Element VonMisses stress.
        RRVector elementVonMisses;

        //! Eigen values.
        RRVector d;
        //! Eigen vectors.
        RRMatrix ev;

        //! Transient analysis is elastodynamic (Newmark / HHT time integration).
        bool dynamic;
        //! Initial acceleration has to be computed from equilibrium.
        bool computeInitialAcceleration;
        //! Newmark gamma.
        double newmarkGamma;
        //! Newmark beta.
        double newmarkBeta;
        //! Hilber-Hughes-Taylor alpha.
        double hhtAlpha;
        //! Time step size.
        double timeStepSize;
        //! Acceleration part of the Newmark predictor on free directions.
        RSparseMatrix D;
        //! Initial equilibrium residual on free directions.
        RRVector r;

    public:

        //! Constructor.
        explicit RSolverStress(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData, bool modalAnalysis);

        //! Destructor.
        ~RSolverStress();

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Number of modes the last eigen-value solve actually produced.
        uint getNComputedModes() const override;

    protected:

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

        //! Solve stress-strain problem.
        void solveStressStrain();

        //! Solve eigen-value problem.
        void solveEigenValue();

        //! Return true if transient analysis is elastodynamic.
        bool findDynamic() const;

        //! Find Newmark and HHT time integration parameters.
        void findTimeIntegrationParameters();

        //! Apply initial displacement and velocity.
        void applyInitialConditions();

        //! Remember the state at the beginning of the time step.
        void storeOldState();

        //! Solve M*a = f - K*u - C*v for the initial acceleration and
        //! correct the right hand side with it.
        void solveInitialAcceleration();

        //! Update velocity and acceleration from the new displacement.
        void updateVelocityAndAcceleration();

        //! Set displacemen.
        void setDisplacement(const RRVector &v);

        //! Convert solution vector in node frames to global node vectors.
        void setNodeVector(const RRVector &v, RSolverCartesianVector<RRVector> &nodeVector, bool prescribed) const;

        //! Store node vector variable.
        void storeNodeVector(RVariableType variableType, const RSolverCartesianVector<RRVector> &nodeVector);

        //! Lump mass matrix - diagonal scaling which preserves element mass.
        static void lumpMassMatrix(RRMatrix &Me);

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Local rotations are built from the constraints themselves, see
        //! generateLocalConstraints(), so the generic geometric pass is not used.
        void updateLocalRotations() override;

        //! Collect every displacement constraint acting on each node, reduce
        //! them to an orthonormal set, and from that build the node local frame,
        //! the number of constrained directions and the prescribed values.
        void generateLocalConstraints();

        //! Generate node book.
        void generateNodeBook();

        //! Return true if given component of the boundary condition is present
        //! and switched on.
        static bool isComponentEnabled(const RBoundaryCondition &bc, RVariableType variableType);

        //! Return the value of given component of the boundary condition at the
        //! current solver time, or zero when the component is absent.
        double findComponentValue(const RBoundaryCondition &bc, RVariableType variableType) const;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, AssemblyBuffer &buffer);

        //! Apply local rotations to matrix.
        void applyLocalRotations(unsigned int elementID, RRMatrix &Ae);

        //! Apply local rotations to vector.
        void applyLocalRotations(unsigned int elementID, RRVector &fe);

};

#endif // RSOLVERSTRESS_H

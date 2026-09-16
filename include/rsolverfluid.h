#ifndef RSOLVERFLUID_H
#define RSOLVERFLUID_H


#include <rbl_stop_watch.h>

#include <rml_element_shape_derivation.h>

#include "rmatrixmanager.h"
#include "rsolvergeneric.h"

class FluidMatrixContainer;

class RSolverFluid : public RSolverGeneric
{

    protected:

        //! How far the residual has to come down within one solve before the
        //! iteration counts as converged, relative to where that solve started.
        //! Guards against a nearly singular system, whose increments are small
        //! because it can not move rather than because it has arrived.
        static const double residualDropRatio;

        //! Step of the central difference used by verifyJacobian().
        static const double differenceStep;

        //! Smallest share of the computed increment the under-relaxation may
        //! fall back to. A step shorter than this buys nothing.
        static const double minRelaxation;
        //! Factor the relaxation is cut by after a pass which made the residual
        //! worse, and recovered by - its reciprocal is not used, recovery is
        //! deliberately slower than retreat - after a pass which improved it.
        static const double relaxationCutFactor;
        static const double relaxationGrowFactor;
        //! How much the residual may rise from one pass to the next before the
        //! rise counts as an overshoot rather than as the wander of a
        //! stabilised formulation.
        static const double relaxationRiseTolerance;

        //! Node height.
        RRVector freePressureNodeHeight;
        //! Element scales.
        RRVector elementScales;

        //! Element pressure.
        RRVector elementPressure;
        //! Element velocity.
        RSolverCartesianVector<RRVector> elementVelocity;
        //! Element acceleration.
        RSolverCartesianVector<RRVector> elementGravity;

        //! Node pressure.
        RRVector nodePressure;
        //! Node velocity.
        RSolverCartesianVector<RRVector> nodeVelocity;
        //! Node velocity.
        RSolverCartesianVector<RRVector> nodeVelocityOld;
        //! Node acceleration.
        RSolverCartesianVector<RRVector> nodeAcceleration;

        //! Stream velocity.
        double streamVelocity;
        double invStreamVelocity;

        //! Element density.
        RRVector elementDensity;
        //! Element viscosity.
        RRVector elementViscosity;

        //! Average density.
        double avgRo;
        //! Average dynamic viscosity.
        double avgU;

        //! Velocity convergence.
        double cvgV;
        //! Pressure convergence.
        double cvgP;

        //! Vector of surface normals.
        std::vector<RR3Vector> elementNormals;
        //! Element gravity magnitude.
        RRVector elementGravityMagnitude;
        //! Vector of element level shape function derivatives.
        std::vector<RElementShapeDerivation *> shapeDerivations;
        //! Cached active local DOF indexes for each assembled element.
        std::vector<std::vector<uint>> elementActiveDofs;
        //! Cached global matrix positions for each element-local matrix entry.
        std::vector<std::vector<uint>> elementMatrixPositions;
        //! Cached global RHS positions for each element-local vector entry.
        std::vector<std::vector<uint>> elementVectorPositions;

        //! Per-thread assembly buffers.
        //! Rebuilt only when the sparse matrix pattern changes.
        std::vector<RSparseMatrix> threadAssemblyMatrices;
        std::vector<RRVector> threadAssemblyVectors;

        //! Stop-watches
        RStopWatch recoveryStopWatch;
        RStopWatch buildStopWatch;
        RStopWatch solverStopWatch;
        RStopWatch updateStopWatch;

        //! Statistics state (replaces function-local statics).
        uint statsCounter;
        double statsOldResidual;
        //! Residual of the current pass.
        double residual;
        //! Residual of the first pass of the current solve, which the residual
        //! of every following pass is measured against.
        double residualFirst;
        //! Residual of the previous pass of the current solve.
        double previousResidual;
        //! Share of the computed increment which is actually applied. One means
        //! the full Newton step.
        double relaxation;

        //! Guards the one-time warm-start of x on a restarted run.
        bool xInitialized;

        //! Whether the finite-difference Jacobian check was asked for.
        static bool verifyJacobianRequested;
        //! Whether it has already run - it runs once.
        bool jacobianVerified;
        //! Whether the element routines take their stabilisation parameters
        //! from the vectors below instead of from the field being assembled.
        bool freezeStabilization;
        //! Element velocity magnitude and direction held fixed while the
        //! Jacobian check perturbs the field. Without them the difference would
        //! also measure the derivative of the stabilisation parameters, which
        //! the assembled matrix deliberately does not contain.
        RRVector frozenMvh;
        std::vector<RR3Vector> frozenS;

    public:

        //! Constructor.
        explicit RSolverFluid(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverFluid() override;

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Ask for the assembled matrix to be compared against a finite
        //! difference of the residual, once, on the next solve.
        //! The check reassembles the whole system twice per degree of freedom,
        //! so it belongs on a mesh of a few elements and nowhere else.
        static void setVerifyJacobian(bool verifyJacobian);

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

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Find input vectors.
        void findInputVectors();

        //! Generate node book.
        void generateNodeBook();

        //! Compute free pressure node height.
        void computeFreePressureNodeHeight();

        //! Compute element shape derivatives.
        void computeShapeDerivatives();

        //! Clear element shape derivatives.
        void clearShapeDerivatives();

        //! Compute element matrix.
        void computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager);

        //! Compute element matrix.
        void computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager);

        //! Compute tetrahedra element matrix.
        void computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager);

        //! Update the residual of the current pass and, from how it compares
        //! with the previous one, the relaxation of the step about to be taken.
        void updateResidualAndRelaxation();

        //! Recompute the nodal acceleration from the current and the previous
        //! velocity field. solve() does this as part of the update; the
        //! Jacobian check has to do it itself, because the residual depends on
        //! the acceleration and the acceleration depends on the velocity being
        //! perturbed.
        void updateNodeAcceleration();

        //! Compare the assembled matrix against a finite difference of the
        //! residual and report where the two disagree.
        void verifyJacobian();

        //! Find time scale.
        double findTimeScale() const;

        //! Find Re scale.
        double findReScale() const;

        //! Find weight scale.
        double findWeightScale() const;

        //! Compute element scales.
        void computeElementScales();

        //! Find element free pressure.
        void computeElementFreePressure(RRVector &values, RBVector &setValues);

        //! Build reusable sparse matrix pattern.
        void buildSparseMatrixPattern(const RBVector &elementFreePressureSetValues);

        //! Assembly matrix.
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Apply local rotations.
        void applyLocalRotations(unsigned int elementID, RRMatrix &Ae);

    public:

        //! Find stream velocity.
        static double computeStreamVelocity(const RModel &rModel,
                                            const RSolverCartesianVector<RRVector> &nodeVelocity,
                                            bool averageBased);

};

#endif // RSOLVERFLUID_H

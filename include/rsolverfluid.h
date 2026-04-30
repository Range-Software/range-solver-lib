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

        //! Stop-watches
        RStopWatch recoveryStopWatch;
        RStopWatch buildStopWatch;
        RStopWatch solverStopWatch;
        RStopWatch updateStopWatch;

        //! Statistics state (replaces function-local statics).
        uint statsCounter;
        double statsOldResidual;

        //! Guards the one-time warm-start of x on a restarted run.
        bool xInitialized;

    public:

        //! Constructor.
        explicit RSolverFluid(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverFluid() override;

        //! Check if solver has converged.
        bool hasConverged() const override;

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

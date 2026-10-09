#ifndef RSOLVERFLUIDPARTICLE_H
#define RSOLVERFLUIDPARTICLE_H

#include <rbl_stop_watch.h>
#include <rml_element_shape_derivation.h>
#include <rml_element_shape_function.h>

#include "rmatrixmanager.h"
#include "rsolvergeneric.h"

class FluidParticleMatrixContainer;

class RSolverFluidParticle : public RSolverGeneric
{

    protected:

        //! Element particle concentration.
        RRVector elementConcentration;
        //! Element particle rate.
        RRVector elementRate;
        //! Element velocity.
        RSolverCartesianVector<RRVector> elementVelocity;

        //! Stream velocity.
        double streamVelocity;

        //! Node particle concentration.
        RRVector nodeConcentration;
        //! Node particle rate.
        RRVector nodeRate;
        //! Node velocity.
        RSolverCartesianVector<RRVector> nodeVelocity;

        //! Element density.
        RRVector elementDensity;
        //! Element density.
        RRVector elementDiffusion;

        //! Maximum saturation (scaled), zero if saturation is not limited.
        double maximumSaturation;

        //! Surface elements on the boundary of the computed volume (a face of exactly one volume element).
        RBVector boundaryFaces;

        //! Concentration convergence.
        double cvgC;

        //! Vector of element level shape function derivatives.
        std::vector<RElementShapeDerivation *> shapeDerivations;

        //! Stop-watches
        RStopWatch recoveryStopWatch;
        RStopWatch buildStopWatch;
        RStopWatch assemblyStopWatch;
        RStopWatch solverStopWatch;
        RStopWatch updateStopWatch;

        //! Statistics state (replaces function-local statics).
        uint statsCounter;
        double statsOldResidual;

    public:

        //! Constructor.
        explicit RSolverFluidParticle(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverFluidParticle() override;

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        //! Generate node rate input vector.
        void generateNodeRateVector();

        //! Return saturation reaction coefficient for given node.
        //! Positive particle rate f is limited by saturation as f * (1 - c / cSat).
        double computeNodeSaturationCoefficient(uint nodeID) const;

        //! Warn if concentration exceeds maximum saturation.
        void checkMaximumSaturation() const;

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

        //! Compute element shape derivatives.
        void computeShapeDerivatives();

        //! Clear element shape derivatives.
        void clearShapeDerivatives();

        //! Find surface elements on the boundary of the computed volume.
        void findBoundaryFaces();

        //! Compute boundary element matrix of flow entering through a free boundary.
        void computeBackflowElement(unsigned int elementID, RRMatrix &Ae, RRVector &be);

        //! Compute element matrix.
        void computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager);

        //! Compute element matrix.
        void computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager);

        //! Compute tetrahedra element matrix.
        void computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager);

        //! Assembly matrix.
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be);

        //! Assembly matrix.
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be, RSparseMatrix &Ap, RRVector &bp);

};

#endif // RSOLVERFLUIDPARTICLE_H

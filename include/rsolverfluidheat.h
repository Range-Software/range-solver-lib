#ifndef RSOLVERFLUIDHEAT_H
#define RSOLVERFLUIDHEAT_H

#include <rbl_stop_watch.h>
#include <rml_element_shape_derivation.h>

#include "rmatrixmanager.h"
#include "rsolvergeneric.h"

class FluidHeatMatrixContainer;

class RSolverFluidHeat : public RSolverGeneric
{

    protected:

        //! Element heat capacity vector.
        RRVector elementCapacity;
        //! Element thermal conduction vector.
        RRVector elementConduction;
        //! Element density vector.
        RRVector elementDensity;
        //! Node temperature.
        RRVector nodeTemperature;
        //! Node heat vector.
        RRVector nodeHeat;
        //! Element heat vector.
        RRVector elementHeat;
        //! Element radiation heat vector.
        RRVector elementRadiativeHeat;
        //! Element joule heat.
        RRVector elementJouleHeat;
        //! Element heat flux vector.
        std::vector<RR3Vector> elementHeatFlux;

        //! Element velocity.
        RSolverCartesianVector<RRVector> elementVelocity;
        //! Node velocity.
        RSolverCartesianVector<RRVector> nodeVelocity;

        //! Stream velocity.
        double streamVelocity;

        //! Temperature convergence.
        double cvgT;

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
        explicit RSolverFluidHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverFluidHeat() override;

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        //! Find temperature scale.
        double findTemperatureScale() const;

        //! Generate node heat input vector.
        void generateNodeHeatVector();

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

        //! Compute element matrix.
        void computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager);

        //! Compute element matrix.
        void computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager);

        //! Compute tetrahedra element matrix.
        void computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager);

        //! Assembly matrix.
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be);

        //! Assembly matrix.
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be, RSparseMatrix &Ap, RRVector &bp);

};

#endif // RSOLVERFLUIDHEAT_H

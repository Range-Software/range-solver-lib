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

        //! Temperature convergence - relative size of the last temperature change.
        double cvgT;

        //! Fluid volume element behind every surface element carrying the Forced
        //! convection condition, or RConstants::eod where there is none.
        RUVector wallFluidElements;
        //! Nodes of the walls - surface elements with a fluid element behind them.
        RBVector wallNodes;
        //! Solid node temperature recovered from the heat solver.
        //! Empty when no heat solve has run.
        RRVector solidNodeTemperature;
        //! Whether the wall nodes are held at the solid temperature.
        bool wallCoupled;
        //! Wall nodes held at the solid temperature in the current pass.
        RBVector coupledWallNodes;
        //! Node temperature the current pass started from - the previous time
        //! level of a transient solve.
        RRVector nodeTemperatureOld;
        //! Wall node temperature the last pass held the walls at - the solid
        //! temperature relaxed. Empty until the first coupled pass.
        RRVector wallTemperature;
        //! Difference between the solid temperature and the wall temperature
        //! of the previous coupled pass. Empty until the second one.
        RRVector wallResidual;
        //! Aitken relaxation factor of the wall temperature.
        double wallRelaxation;
        //! Wall heat transfer coefficient - the conductance of the first fluid
        //! element. Negative on every element which is not a wall.
        RRVector elementWallHtc;
        //! Wall reference temperature - the fluid temperature the conductance acts
        //! across. Meaningful only where elementWallHtc is not negative.
        RRVector elementWallHtt;

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

        //! Key the wall heat transfer coefficient is shared under, so the heat
        //! solver can drive its Forced convection walls with it. An element
        //! vector, negative on every element which is not a wall.
        static const QString wallHeatTransferCoefficientKey;

        //! Key the wall reference temperature is shared under, for the same
        //! reason. An element vector paired with the coefficient.
        static const QString wallFluidTemperatureKey;

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

        //! Find the fluid volume element behind every surface element carrying
        //! the Forced convection condition, and the nodes of such walls.
        void findWallElements();

        //! Hold the wall nodes at the temperature the heat solver computed in
        //! the solid, once it has run. Until then the walls are adiabatic.
        //! The temperature is relaxed with the Aitken factor, which takes the
        //! alternation of the two solves to the coupled solution in a few passes
        //! where plain alternation can take hundreds.
        void applyWallTemperature();

        //! Compute the wall heat transfer coefficient and reference temperature.
        //! Together they reproduce the heat flux the fluid solve takes through
        //! the wall.
        void computeWallHeatTransfer();

        //! Return the heat entering the fluid through every wall node held at
        //! the solid temperature - the residual of the fluid system at that node,
        //! which is the flux consistent with the discretisation.
        RRVector computeWallReaction();

        //! Store solver results into the shared data container.
        void storeSharedData() override;

        //! Recover previously computed results from the shared data container.
        void recoverSharedData() override;

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

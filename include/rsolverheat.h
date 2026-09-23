#ifndef RSOLVERHEAT_H
#define RSOLVERHEAT_H

#include <rbl_r3vector.h>

#include "rsolvergeneric.h"

class RSolverHeat : public RSolverGeneric
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
        //! Element heat vector.
        RRVector elementHeat;
        //! Element heat rate per unit area - surface elements only.
        RRVector elementHeatRateArea;
        //! Element heat rate per unit volume - volume elements only.
        RRVector elementHeatRateVolume;
        //! Element radiation heat vector.
        RRVector elementRadiativeHeat;
        //! Element joule heat.
        RRVector elementJouleHeat;
        //! Element heat flux vector.
        std::vector<RR3Vector> elementHeatFlux;
        //! Wall heat transfer coefficient recovered from the fluid heat solver,
        //! negative on every element which is not a wall.
        //! Empty when no fluid heat solve has run.
        RRVector fluidWallHtc;
        //! Wall reference temperature recovered from the fluid heat solver.
        //! Empty when no fluid heat solve has run.
        RRVector fluidWallHtt;
        //! Whether any wall took its convection from the fluid heat solver.
        bool wallCoupled;
        //! Temperature convergence - relative size of the last temperature change.
        double cvgT;
        //! Element heat transfer coefficient.
        RRVector elementHeatTransferCoefficient;

    public:

        //! Key the solved node temperature is shared under, so the fluid heat
        //! solver can hold its walls at the temperature of the solid.
        static const QString solidNodeTemperatureKey;

        //! Constructor.
        explicit RSolverHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverHeat();

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        //! Store solver results into the shared data container.
        void storeSharedData() override;

        //! Recover previously computed results from the shared data container.
        void recoverSharedData() override;

        //! Heat conduction is solved in solids only. A fluid domain is left to the
        //! fluid heat solver, and so is any boundary entity lying inside one.
        void findComputableElements(RProblemType problemType) override;

        //! Return the heat transfer coefficient and the reference temperature the
        //! fluid heat solver computed for the given wall element. False when no
        //! such result covers it - the surface borders no fluid, or no fluid heat
        //! solve has run yet.
        bool findFluidWall(unsigned int elementId, double &htc, double &htt) const;

        //! Find temperature scale.
        double findTemperatureScale() const;

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

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Get simple convection BC values.
        bool getSimpleConvection(const RElementGroup &elementGroup, double &htc, double &htt);

        //! Throw when a value configured on a convection condition would leave the
        //! correlation with nothing to work with - the dimensionless groups divide by
        //! the viscosity, the thermal conductivity and the hydraulic diameter, and a
        //! zero density, heat capacity or mean velocity collapses them just as surely.
        void checkConvectionInput(double value,
                                  RVariableType variableType,
                                  RBoundaryConditionType boundaryConditionType,
                                  const RElementGroup &elementGroup) const;

        //! Report where a Forced convection condition takes its heat transfer
        //! from - the fluid heat solver, the correlation, or neither.
        void reportForcedConvection(const RElementGroup &elementGroup);

        //! Get forced convection BC values. A wall bordering a fluid domain takes
        //! them from the fluid heat solver results. Only where those are missing
        //! are they correlated from the values configured on the condition.
        bool getForcedConvection(const RElementGroup &elementGroup, unsigned int elementId, double &htc, double &htt);

        //! Get natural convection BC values.
        bool getNaturalConvection(const RElementGroup &elementGroup, unsigned int elementId, double &htc, double &htt);

};

#endif // RSOLVERHEAT_H

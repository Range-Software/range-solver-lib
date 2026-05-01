#ifndef RSOLVERELECTROSTATICS_H
#define RSOLVERELECTROSTATICS_H

#include "rsolvergeneric.h"

class RSolverElectrostatics : public RSolverGeneric
{

    protected:

        //! Node electric potential.
        RRVector nodeElectricPotential;
        //! Element electric field.
        std::vector<RR3Vector> elementElectricField;
        //! Element current density.
        std::vector<RR3Vector> elementCurrentDensity;
        //! Element electric energy.
        RRVector elementElectricEnergy;
        //! Element electric resistivity.
        RRVector elementElectricResistivity;
        //! Element joule heat.
        RRVector elementJouleHeat;
        //! Element relative permittivity.
        RRVector elementRelativePermittivity;
        //! Element electric conductivity.
        RRVector elementElectricConductivity;

    public:

        //! Constructor.
        explicit RSolverElectrostatics(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverElectrostatics();

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

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ke, const RRVector &fe);

};

#endif // RSOLVERELECTROSTATICS_H

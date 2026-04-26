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
        bool hasConverged() const;

    protected:

        //! Update scales.
        void updateScales();

        //! Recover previously computed results.
        void recover();

        //! Prepare solver.
        void prepare();

        //! Run matrix solver.
        void solve();

        //! Process solver results.
        void process();

        //! Store solver results.
        void store();

        //! Process statistics.
        void statistics();

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ke, const RRVector &fe);

};

#endif // RSOLVERELECTROSTATICS_H

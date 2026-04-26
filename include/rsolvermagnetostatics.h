#ifndef RSOLVERMAGNETOSTATICS_H
#define RSOLVERMAGNETOSTATICS_H

#include "rsolvergeneric.h"

class RSolverMagnetostatics : public RSolverGeneric
{

    protected:

        //! Node current density.
        RSolverCartesianVector<RRVector> nodeCurrentDensity;
        //! Node current density.
        RSolverCartesianVector<RRVector> nodeMagneticField;

    public:

        //! Constructor.
        explicit RSolverMagnetostatics(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverMagnetostatics();

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

#endif // RSOLVERMAGNETOSTATICS_H

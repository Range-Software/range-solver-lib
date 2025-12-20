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
        bool hasConverged(void) const;

    protected:

        //! Update scales.
        void updateScales(void);

        //! Recover previously computed results.
        void recover(void);

        //! Prepare solver.
        void prepare(void);

        //! Run matrix solver.
        void solve(void);

        //! Process solver results.
        void process(void);

        //! Store solver results.
        void store(void);

        //! Process statistics.
        void statistics(void);

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ke, const RRVector &fe);

};

#endif // RSOLVERMAGNETOSTATICS_H

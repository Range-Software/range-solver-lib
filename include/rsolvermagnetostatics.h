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

#endif // RSOLVERMAGNETOSTATICS_H

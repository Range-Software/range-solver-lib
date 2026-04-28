#ifndef RSOLVERWAVE_H
#define RSOLVERWAVE_H

#include "rsolvergeneric.h"

class RSolverWave : public RSolverGeneric
{

    protected:

        //! Element wave speed.
        RRVector elementWaveSpeed;
        //! Element wave displacement.
        RRVector elementWaveDisplacement;
        //! Element node displacement.
        RRVector nodeWaveDisplacement;

    public:

        //! Constructor.
        explicit RSolverWave(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverWave();

        //! Check if solver has converged.
        bool hasConverged() const;

    protected:

        //! Initialize solver.
        void initialize() override;

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

};

#endif // RSOLVERWAVE_H

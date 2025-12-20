#ifndef RSOLVERMESH_H
#define RSOLVERMESH_H

#include "rsolvergeneric.h"

class RSolverMesh : public RSolverGeneric
{

    protected:

        //! Mesh input.
        RMeshInput meshInput;

    public:

        //! Constructor.
        explicit RSolverMesh(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverMesh();

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

};

#endif // RSOLVERMESH_H

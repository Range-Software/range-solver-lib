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

#endif // RSOLVERMESH_H

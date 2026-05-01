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

};

#endif // RSOLVERMESH_H

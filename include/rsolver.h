#ifndef RSOLVER_H
#define RSOLVER_H

#include <rml_model.h>

#include "rsolvergeneric.h"
#include "rsolvershareddata.h"

class RSolver
{

    protected:

        //! Pointer to model object.
        RModel *pModel;
        //! Model file name.
        QString modelFileName;
        //! Convergence file name.
        QString convergenceFileName;
        //! Shared variables.
        RSolverSharedData sharedData;
        //! Map of solvers.
        QMap<RProblemTypeMask,RSolverGeneric*> solvers;
        //! Map of sover type and execution count.
        QMap<RProblemType,uint> solversExecutionCount;

    private:

        //! Internal initialization function.
        void _init(const RSolver *pSolver = nullptr);

    public:

        //! Constructor.
        RSolver(RModel &model, const QString &modelFileName, const QString &convergenceFileName);

        //! Copy constructor.
        RSolver(const RSolver &solver);

        //! Destructor.
        ~RSolver();

        //! Assignment operator.
        RSolver & operator =(const RSolver &solver);

        //! Run solver.
        void run(void);

    protected:

        //! Run single solver.
        void runSingle(void);

        //! Run single solver for given problem task.
        //! Return task convergence status (true = converged).
        bool runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration);

};

#endif // RSOLVER_H

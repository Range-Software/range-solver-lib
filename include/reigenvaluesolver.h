#ifndef REIGENVALUESOLVER_H
#define REIGENVALUESOLVER_H

#include <rbl_utils.h>
#include <rml_eigen_value_solver_conf.h>

#include "rmatrixsolver.h"

class REigenValueSolver
{

    protected:

        //! Matrix solver input.
        REigenValueSolverConf eigenValueSolverConf;

        //! Matrix solver input.
        RMatrixSolverConf matrixSolverConf;

    private:

        //! Internal initialization function.
        void _init(const REigenValueSolver *pEigenValueSolver = nullptr);

    public:

        //! Constructor.
        REigenValueSolver(const REigenValueSolverConf &eigenValueSolverConf, const RMatrixSolverConf &matrixSolverConf);

        //! Copy constructor.
        REigenValueSolver(const REigenValueSolver &eigenValueSolver);

        //! Destructor.
        ~REigenValueSolver();

        //! Assignment operator.
        REigenValueSolver &operator =(const REigenValueSolver &eigenValueSolver);

        //! Solve matrix system.
        //! d = eigen values
        //! ev = eigen vectors
        void solve(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev);

    protected:

        //! Lanczos method solver.
        void solveLanczos(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d);

        //! Arnoldi method.
        void solveArnoldi(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev);

        //! Rayleigh method.
        void solveRayleigh(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &e, RRMatrix &ev);

        //! QL decomposition.
        static void qlDecomposition(RRVector &d, RRVector &e);

        //! QR decomposition - Gram-Schmidt process
        static void qrDecomposition(const RRMatrix &H, uint nIterations, double convergenceValue, RRVector &d, RRMatrix &V);

};

#endif // REIGENVALUESOLVER_H

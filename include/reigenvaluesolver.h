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

        //! Subspace iteration with a Rayleigh-Ritz projection.
        //! Converges towards the lowest eigen values of K*phi = lambda*M*phi.
        void solveSubspaceIteration(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev);

        //! Inverse power iteration with a Rayleigh quotient.
        //! Converges towards the single lowest eigen pair.
        void solveInversePowerIteration(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev);

        //! QL decomposition.
        static void qlDecomposition(RRVector &d, RRVector &e);

        //! Eigen values and vectors of a small dense symmetric matrix, using
        //! the cyclic Jacobi method. A is destroyed, the eigen vectors are the
        //! columns of V.
        static void jacobiEigen(RRMatrix &A, RRVector &d, RRMatrix &V);

        //! Cholesky decomposition of a small dense symmetric positive definite
        //! matrix, A = L * L^T. Returns false if A is not positive definite.
        static bool choleskyDecomposition(const RRMatrix &A, RRMatrix &L);

        //! Solve L*x = b for a lower triangular L.
        static void forwardSubstitution(const RRMatrix &L, const RRVector &b, RRVector &x);

        //! Solve L^T*x = b for a lower triangular L.
        static void backwardSubstitution(const RRMatrix &L, const RRVector &b, RRVector &x);

        //! Orthonormalize the columns of X with a modified Gram-Schmidt pass.
        //! Returns the number of independent columns kept.
        static uint orthonormalizeColumns(RRMatrix &X);

};

#endif // REIGENVALUESOLVER_H

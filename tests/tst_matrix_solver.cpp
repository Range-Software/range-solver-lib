#include <QtTest>

#include "rmatrixsolver.h"
#include "rml_matrix_solver_conf.h"
#include "rml_sparse_matrix.h"
#include "rbl_rvector.h"

namespace
{
constexpr double tol = 1.0e-6;

RMatrixSolverConf makeCgConf()
{
    RMatrixSolverConf conf(RMatrixSolverConf::CG);
    conf.setNInnerIterations(1000);
    conf.setNOuterIterations(1000);
    conf.setSolverCvgValue(1.0e-13);
    conf.setOutputFrequency(0);
    return conf;
}
}

class TestMatrixSolver : public QObject
{
    Q_OBJECT

private slots:

    void solveSymmetricPositiveDefinite();
    void solveDiagonalSystem();
};

void TestMatrixSolver::solveSymmetricPositiveDefinite()
{
    // A = [[4,1],[1,3]], b = [1,2]  ->  x = [1/11, 7/11].
    RSparseMatrix a;
    a.setNRows(2);
    a.addValue(0, 0, 4.0);
    a.addValue(0, 1, 1.0);
    a.addValue(1, 0, 1.0);
    a.addValue(1, 1, 3.0);

    RRVector b(2, 0.0);
    b[0] = 1.0;
    b[1] = 2.0;

    RRVector x;
    RMatrixSolver solver(makeCgConf());
    solver.solve(a, b, x);

    QCOMPARE(x.getNRows(), 2u);
    QVERIFY(qAbs(x[0] - 1.0 / 11.0) < tol);
    QVERIFY(qAbs(x[1] - 7.0 / 11.0) < tol);
}

void TestMatrixSolver::solveDiagonalSystem()
{
    // Diagonal A trivially gives x[i] = b[i] / A[i][i].
    RSparseMatrix a;
    a.setNRows(3);
    a.addValue(0, 0, 2.0);
    a.addValue(1, 1, 4.0);
    a.addValue(2, 2, 5.0);

    RRVector b(3, 0.0);
    b[0] = 6.0;
    b[1] = 8.0;
    b[2] = 10.0;

    RRVector x;
    RMatrixSolver solver(makeCgConf());
    solver.solve(a, b, x);

    QVERIFY(qAbs(x[0] - 3.0) < tol);
    QVERIFY(qAbs(x[1] - 2.0) < tol);
    QVERIFY(qAbs(x[2] - 2.0) < tol);
}

QTEST_APPLESS_MAIN(TestMatrixSolver)

#include "tst_matrix_solver.moc"

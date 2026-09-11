#include <QtTest>

#include <cmath>

#include "reigenvaluesolver.h"
#include "rml_matrix_solver_conf.h"
#include "rml_sparse_matrix.h"

namespace
{

RMatrixSolverConf makeCgConf()
{
    RMatrixSolverConf conf(RMatrixSolverConf::CG);
    conf.setNInnerIterations(200);
    conf.setNOuterIterations(2000);
    conf.setSolverCvgValue(1.0e-14);
    conf.setOutputFrequency(0);
    return conf;
}

REigenValueSolverConf makeEigenConf(REigenValueSolverConf::Method method, uint nEigenValues)
{
    REigenValueSolverConf conf(method);
    conf.setNEigenValues(nEigenValues);
    conf.setNIterations(500);
    conf.setSolverCvgValue(1.0e-14);
    conf.setOutputFrequency(0);
    return conf;
}

//! K = [[2,-1],[-1,1]] * 1e6, M = identity.
//!
//! The generalised eigenvalues lambda are around 3.8e5 and 2.6e6, so lambda and
//! its reciprocal are more than ten orders of magnitude apart. That makes it
//! unambiguous which of the two a returned value is on, without relying on the
//! accuracy of the iteration itself.
void buildSystem(RSparseMatrix &m, RSparseMatrix &k)
{
    const double scale = 1.0e6;

    m.clear();
    m.setNRows(2);
    m.addValue(0,0,1.0);
    m.addValue(1,1,1.0);

    k.clear();
    k.setNRows(2);
    k.addValue(0,0,2.0*scale);
    k.addValue(0,1,-1.0*scale);
    k.addValue(1,0,-1.0*scale);
    k.addValue(1,1,1.0*scale);
}

}

class TestEigenValueSolver : public QObject
{
    Q_OBJECT

private slots:

    void multipleValuesAreGeneralisedEigenValues();
    void singleValueUsesTheSameScale();
    void inversePowerIterationFindsTheLowestValue();
    void valuesAreSortedAscending();
};

//! Several extracted values have to come back as lambda of
//! K * phi = lambda * M * phi, not as its reciprocal.
void TestEigenValueSolver::multipleValuesAreGeneralisedEigenValues()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,2),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),2u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    double lambdaHigh = (3.0 + std::sqrt(5.0)) / 2.0 * 1.0e6;

    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Lowest eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));
    QVERIFY2(std::fabs(d[1]-lambdaHigh) < 1.0e-6*lambdaHigh,
             qPrintable(QString("Highest eigen value = %1, expected %2").arg(d[1]).arg(lambdaHigh)));
}

//! Extracting a single value used to skip the inversion, handing back a value
//! on the reciprocal scale. One and several values have to agree in scale.
void TestEigenValueSolver::singleValueUsesTheSameScale()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,1),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),1u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Single eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));
}

//! The inverse power iteration has to converge to the same lowest eigen value.
void TestEigenValueSolver::inversePowerIterationFindsTheLowestValue()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::InversePowerIteration,1),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),1u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));

    QCOMPARE(ev.getNRows(),1u);
    QCOMPARE(ev.getNColumns(),2u);
}

//! The lowest eigen value has to come first, so that mode 0 is the fundamental.
void TestEigenValueSolver::valuesAreSortedAscending()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,2),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),2u);
    QVERIFY2(d[0] <= d[1],
             qPrintable(QString("Eigen values %1, %2 are not sorted ascending").arg(d[0]).arg(d[1])));
}

QTEST_APPLESS_MAIN(TestEigenValueSolver)

#include "tst_eigen_value_solver.moc"

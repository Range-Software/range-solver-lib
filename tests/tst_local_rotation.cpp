#include <QtTest>

#include "rlocalrotation.h"
#include "rbl_rmatrix.h"
#include "rbl_rvector.h"

namespace
{
constexpr double tol = 1.0e-9;

// Rotation by +90 degrees about the Z axis.
RRMatrix rotationZ90()
{
    RRMatrix r(3, 3, 0.0);
    r[0][0] = 0.0; r[0][1] = -1.0; r[0][2] = 0.0;
    r[1][0] = 1.0; r[1][1] = 0.0;  r[1][2] = 0.0;
    r[2][0] = 0.0; r[2][1] = 0.0;  r[2][2] = 1.0;
    return r;
}
}

class TestLocalRotation : public QObject
{
    Q_OBJECT

private slots:

    void inactiveByDefault();
    void activateDeactivate();
    void inverseIsComputed();
    void rotateVector();
};

void TestLocalRotation::inactiveByDefault()
{
    RLocalRotation rotation;
    QVERIFY(!rotation.isActive());
}

void TestLocalRotation::activateDeactivate()
{
    RLocalRotation rotation;
    rotation.activate(rotationZ90());
    QVERIFY(rotation.isActive());
    rotation.deactivate();
    QVERIFY(!rotation.isActive());
}

void TestLocalRotation::inverseIsComputed()
{
    RLocalRotation rotation;
    rotation.activate(rotationZ90());

    // R * iR must be the identity.
    RRMatrix product;
    RRMatrix::mlt(rotation.getR(), rotation.getInverseR(), product);
    for (uint i = 0; i < 3; i++)
    {
        for (uint j = 0; j < 3; j++)
        {
            QVERIFY(qAbs(product.getValue(i, j) - (i == j ? 1.0 : 0.0)) < tol);
        }
    }
}

void TestLocalRotation::rotateVector()
{
    RLocalRotation rotation;
    rotation.activate(rotationZ90());

    RRVector v(1.0, 0.0, 0.0);
    rotation.rotateResultsVector(v);

    // (1,0,0) rotated by +90 deg about Z -> (0,1,0).
    QVERIFY(qAbs(v[0]) < tol);
    QVERIFY(qAbs(v[1] - 1.0) < tol);
    QVERIFY(qAbs(v[2]) < tol);
}

QTEST_APPLESS_MAIN(TestLocalRotation)

#include "tst_local_rotation.moc"

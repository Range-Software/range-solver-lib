#include <QtTest>

#include "rscales.h"

class TestScales : public QObject
{
    Q_OBJECT

private slots:

    void defaultsAreUnity();
    void settersAndGetters();
    void copyAndAssignment();
};

void TestScales::defaultsAreUnity()
{
    RScales scales;
    QCOMPARE(scales.getMetre(), 1.0);
    QCOMPARE(scales.getKilogram(), 1.0);
    QCOMPARE(scales.getSecond(), 1.0);
    QCOMPARE(scales.getAmpere(), 1.0);
    QCOMPARE(scales.getKelvin(), 1.0);
    QCOMPARE(scales.getCandela(), 1.0);
    QCOMPARE(scales.getMole(), 1.0);
}

void TestScales::settersAndGetters()
{
    RScales scales;
    scales.setMetre(2.0);
    scales.setKilogram(3.0);
    scales.setSecond(4.0);
    scales.setKelvin(5.0);

    QCOMPARE(scales.getMetre(), 2.0);
    QCOMPARE(scales.getKilogram(), 3.0);
    QCOMPARE(scales.getSecond(), 4.0);
    QCOMPARE(scales.getKelvin(), 5.0);
}

void TestScales::copyAndAssignment()
{
    RScales scales;
    scales.setMetre(2.5);
    scales.setMole(6.0);

    RScales copy(scales);
    QCOMPARE(copy.getMetre(), 2.5);
    QCOMPARE(copy.getMole(), 6.0);

    RScales assigned;
    assigned = scales;
    QCOMPARE(assigned.getMetre(), 2.5);
    QCOMPARE(assigned.getMole(), 6.0);

    // Modifying the original must not affect the copy.
    scales.setMetre(9.0);
    QCOMPARE(copy.getMetre(), 2.5);
}

QTEST_APPLESS_MAIN(TestScales)

#include "tst_scales.moc"

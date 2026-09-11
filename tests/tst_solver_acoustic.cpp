#include <QtTest>

#include <cmath>

#include "rml_model.h"
#include "rsolver.h"

namespace
{

// Air.
constexpr double density = 1.2;      // [kg/m^3]
constexpr double soundSpeed = 340.0; // [m/s]

constexpr double ductLength = 1.0;   // [m]
constexpr uint   nElements = 200;

//! Build a one dimensional duct discretised with truss elements.
//!
//! The duct is a single line entity with a unit cross area so that the domain
//! and the boundary integrals share the same measure. A point entity at each
//! end carries the boundary conditions.
RModel buildDuct(double velocityAmplitude, double absorptionCoefficient)
{
    RModel model;

    model.setNNodes(nElements+1);
    for (uint i=0;i<=nElements;i++)
    {
        model.getNode(i).set(ductLength*double(i)/double(nElements),0.0,0.0);
    }

    // Line elements forming the duct. RModel::addElement always files the
    // element into an entity group, so let it build the groups directly -
    // adding a second group with the same elements would assemble them twice.
    for (uint i=0;i<nElements;i++)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,i);
        element.setNodeId(1,i+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Gas);
    material.setName("Air");

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);

    RMaterialProperty soundSpeedProperty(RMaterialProperty::SoundSpeed);
    soundSpeedProperty.add(293.15,soundSpeed);
    material.add(soundSpeedProperty);

    RLine &line = model.getLine(0);
    line.setName("Duct");
    line.setCrossArea(1.0);
    line.setMaterial(material);

    // Driven end.
    RElement inletElement(R_ELEMENT_POINT);
    inletElement.setNodeId(0,0);
    model.addElement(inletElement,true,0);

    RPoint &inlet = model.getPoint(0);
    inlet.setName("Source");
    inlet.setVolume(0.0);

    RBoundaryCondition velocityBc(R_BOUNDARY_CONDITION_VELOCITY_NORMAL);
    velocityBc.getComponent(velocityBc.findComponentPosition(R_VARIABLE_VELOCITY)).add(0.0,velocityAmplitude);
    inlet.addBoundaryCondition(velocityBc);

    // Terminated end.
    RElement outletElement(R_ELEMENT_POINT);
    outletElement.setNodeId(0,nElements);
    model.addElement(outletElement,true,1);

    RPoint &outlet = model.getPoint(1);
    outlet.setName("Termination");
    outlet.setVolume(0.0);

    RBoundaryCondition absorbingBc(R_BOUNDARY_CONDITION_ABSORBING_BOUNDARY);
    absorbingBc.getComponent(absorbingBc.findComponentPosition(R_VARIABLE_ACOUSTIC_ABSORPTION_COEFFICIENT)).add(0.0,absorptionCoefficient);
    outlet.addBoundaryCondition(absorbingBc);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_ACOUSTICS));

    // Keep the stock matrix solver settings - the solver has to cope with them
    // on its own - only silence the per iteration logging.
    model.getMatrixSolverConf(RMatrixSolverConf::GMRES).setOutputFrequency(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    return model;
}

double nodeValue(const RModel &model, RVariableType variableType, uint nodeID)
{
    uint position = model.findVariable(variableType);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(0,nodeID);
}

}

class TestSolverAcoustic : public QObject
{
    Q_OBJECT

private slots:

    void harmonicPlaneWaveInAnechoicDuct();
    void harmonicRigidDuctResonance();
    void transientKeepsNewmarkState();
};

//! A duct driven at one end and anechoically terminated at the other carries a
//! plane travelling wave, for which |p| = rho * c * v everywhere.
void TestSolverAcoustic::harmonicPlaneWaveInAnechoicDuct()
{
    const double velocityAmplitude = 0.01; // [m/s]

    RModel model = buildDuct(velocityAmplitude,1.0);

    RAcousticSetup &acousticSetup = model.getProblemSetup().getAcousticSetup();
    acousticSetup.setAnalysisType(R_ACOUSTIC_ANALYSIS_HARMONIC);
    acousticSetup.setFrequencyStart(500.0);
    acousticSetup.setFrequencyStep(0.0);
    acousticSetup.setNFrequencies(1);

    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    double expected = density * soundSpeed * velocityAmplitude;

    for (uint i=0;i<=nElements;i++)
    {
        double p = nodeValue(model,R_VARIABLE_ACOUSTIC_PRESSURE,i);
        QVERIFY2(std::isfinite(p),qPrintable(QString("Non finite pressure at node %1").arg(i)));
        QVERIFY2(std::fabs(p-expected) < 0.01*expected,
                 qPrintable(QString("Node %1: |p| = %2, expected %3").arg(i).arg(p).arg(expected)));
    }

    // The velocity potential amplitude of the travelling wave is v / k.
    double k = acousticSetup.getAngularFrequency() / soundSpeed;
    double phiExpected = velocityAmplitude / k;
    for (uint i=0;i<=nElements;i+=25)
    {
        double phiRe = nodeValue(model,R_VARIABLE_POTENTIAL,i);
        double phiIm = nodeValue(model,R_VARIABLE_POTENTIAL_IMAGINARY,i);
        double phi = std::sqrt(phiRe*phiRe + phiIm*phiIm);
        QVERIFY2(std::fabs(phi-phiExpected) < 0.01*phiExpected,
                 qPrintable(QString("Node %1: |phi| = %2, expected %3").arg(i).arg(phi).arg(phiExpected)));
    }

    // Sound pressure level of the RMS pressure against the default reference.
    double splExpected = 20.0 * std::log10(expected/std::sqrt(2.0)/R_ACOUSTIC_REFERENCE_PRESSURE);
    double spl = nodeValue(model,R_VARIABLE_ACOUSTIC_SOUND_PRESSURE_LEVEL,nElements/2);
    QVERIFY2(std::fabs(spl-splExpected) < 0.1,
             qPrintable(QString("SPL = %1 dB, expected %2 dB").arg(spl).arg(splExpected)));

    // Time averaged intensity of a plane wave is 0.5 * |p| * |u|.
    uint intensityPosition = model.findVariable(R_VARIABLE_ACOUSTIC_INTENSITY);
    QVERIFY(intensityPosition != RConstants::eod);
    double intensity = model.getVariable(intensityPosition).getValue(0,nElements/2);
    double intensityExpected = 0.5 * expected * velocityAmplitude;
    QVERIFY2(std::fabs(intensity-intensityExpected) < 0.02*intensityExpected,
             qPrintable(QString("Intensity = %1 W/m^2, expected %2 W/m^2").arg(intensity).arg(intensityExpected)));

    // The particle velocity of a plane wave equals the driving velocity.
    uint velocityPosition = model.findVariable(R_VARIABLE_ACOUSTIC_PARTICLE_VELOCITY);
    QVERIFY(velocityPosition != RConstants::eod);
    double particleVelocity = model.getVariable(velocityPosition).getValue(0,nElements/2);
    QVERIFY2(std::fabs(particleVelocity-velocityAmplitude) < 0.01*velocityAmplitude,
             qPrintable(QString("Particle velocity = %1 m/s, expected %2 m/s").arg(particleVelocity).arg(velocityAmplitude)));

    // The phase must advance linearly with the wave number k = omega / c.
    double x1 = ductLength * 50.0 / double(nElements);
    double x2 = ductLength * 150.0 / double(nElements);
    double phase1 = nodeValue(model,R_VARIABLE_ACOUSTIC_PHASE,50);
    double phase2 = nodeValue(model,R_VARIABLE_ACOUSTIC_PHASE,150);

    double dPhase = phase1 - phase2;
    while (dPhase > 180.0) dPhase -= 360.0;
    while (dPhase < -180.0) dPhase += 360.0;

    double dPhaseExpected = k * (x2 - x1) * 180.0 / RConstants::pi;
    while (dPhaseExpected > 180.0) dPhaseExpected -= 360.0;
    while (dPhaseExpected < -180.0) dPhaseExpected += 360.0;

    QVERIFY2(std::fabs(dPhase-dPhaseExpected) < 5.0,
             qPrintable(QString("Phase advance %1 deg, expected %2 deg").arg(dPhase).arg(dPhaseExpected)));
}

//! With a rigid termination the duct resonates. Driving it away from a
//! resonance must still give a bounded, finite standing wave whose amplitude
//! exceeds the travelling wave amplitude.
void TestSolverAcoustic::harmonicRigidDuctResonance()
{
    const double velocityAmplitude = 0.01;

    RModel model = buildDuct(velocityAmplitude,0.0);

    RAcousticSetup &acousticSetup = model.getProblemSetup().getAcousticSetup();
    acousticSetup.setAnalysisType(R_ACOUSTIC_ANALYSIS_HARMONIC);
    // Quarter wave resonance of a 1 m duct is at c/(4L) = 85 Hz - stay away.
    acousticSetup.setFrequencyStart(200.0);
    acousticSetup.setFrequencyStep(0.0);
    acousticSetup.setNFrequencies(1);

    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    double travelling = density * soundSpeed * velocityAmplitude;

    for (uint i=0;i<=nElements;i+=10)
    {
        double p = nodeValue(model,R_VARIABLE_ACOUSTIC_PRESSURE,i);
        QVERIFY2(std::isfinite(p),qPrintable(QString("Non finite pressure at node %1").arg(i)));
    }

    double pMax = 0.0;
    for (uint i=0;i<=nElements;i++)
    {
        pMax = std::max(pMax,nodeValue(model,R_VARIABLE_ACOUSTIC_PRESSURE,i));
    }
    QVERIFY2(pMax > travelling,
             qPrintable(QString("Standing wave peak %1 should exceed travelling wave %2").arg(pMax).arg(travelling)));
}

//! The Newmark state has to survive from one time step to the next, otherwise
//! the transient integration silently degenerates.
void TestSolverAcoustic::transientKeepsNewmarkState()
{
    RModel model = buildDuct(0.01,1.0);

    model.getProblemSetup().getAcousticSetup().setAnalysisType(R_ACOUSTIC_ANALYSIS_TRANSIENT);

    RTimeSolver &timeSolver = model.getTimeSolver();
    timeSolver.setEnabled(true);
    timeSolver.setTimes(RTimeSolver::findTimesVector(20,0.0,1.0e-5));
    timeSolver.setInputNTimeSteps(20);
    timeSolver.setInputStartTime(0.0);
    timeSolver.setInputTimeStepSize(1.0e-5);
    timeSolver.setOutputFrequency(0);

    RSolver solver(model,QString(),QString());
    solver.run();

    // The solver has to publish the Newmark state so that a restart or the next
    // time step can pick it up again.
    QVERIFY(model.findVariable(R_VARIABLE_POTENTIAL_VELOCITY) != RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_POTENTIAL_ACCELERATION) != RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_ACOUSTIC_SOUND_PRESSURE_LEVEL) != RConstants::eod);

    bool nonZeroPotential = false;
    bool allFinite = true;
    for (uint i=0;i<=nElements;i++)
    {
        double phi = nodeValue(model,R_VARIABLE_POTENTIAL,i);
        double p = nodeValue(model,R_VARIABLE_ACOUSTIC_PRESSURE,i);
        allFinite = allFinite && std::isfinite(phi) && std::isfinite(p);
        if (std::fabs(phi) > 1.0e-12)
        {
            nonZeroPotential = true;
        }
    }
    QVERIFY2(allFinite,"Transient solution contains non finite values");
    QVERIFY2(nonZeroPotential,"Velocity source did not excite the duct");

    // A driven, anechoically terminated duct must not blow up.
    double pMax = 0.0;
    for (uint i=0;i<=nElements;i++)
    {
        pMax = std::max(pMax,std::fabs(nodeValue(model,R_VARIABLE_ACOUSTIC_PRESSURE,i)));
    }
    QVERIFY2(pMax < 100.0 * density * soundSpeed * 0.01,
             qPrintable(QString("Transient solution diverged, max |p| = %1").arg(pMax)));
}

QTEST_APPLESS_MAIN(TestSolverAcoustic)

#include "tst_solver_acoustic.moc"

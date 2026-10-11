#include <QtTest>

#include <cmath>

#include "rml_model.h"
#include "rsolver.h"

namespace
{

constexpr double density = 1.2;   // [kg/m^3]
constexpr uint   nDivisions = 3;  // elements along each side of the unit cube

//! Add tetrahedron with positive orientation so that all elements share the same Jacobian sign.
void addTetra(RModel &model, uint groupID, uint n0, uint n1, uint n2, uint n3)
{
    RR3Vector p0(model.getNode(n0).getX(),model.getNode(n0).getY(),model.getNode(n0).getZ());
    RR3Vector p1(model.getNode(n1).getX(),model.getNode(n1).getY(),model.getNode(n1).getZ());
    RR3Vector p2(model.getNode(n2).getX(),model.getNode(n2).getY(),model.getNode(n2).getZ());
    RR3Vector p3(model.getNode(n3).getX(),model.getNode(n3).getY(),model.getNode(n3).getZ());

    RR3Vector a(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]);
    RR3Vector b(p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]);
    RR3Vector c(p3[0]-p0[0],p3[1]-p0[1],p3[2]-p0[2]);

    double volume = a[0]*(b[1]*c[2]-b[2]*c[1])
                  - a[1]*(b[0]*c[2]-b[2]*c[0])
                  + a[2]*(b[0]*c[1]-b[1]*c[0]);
    if (volume < 0.0)
    {
        std::swap(n1,n2);
    }

    RElement element(R_ELEMENT_TETRA1);
    element.setNodeId(0,n0);
    element.setNodeId(1,n1);
    element.setNodeId(2,n2);
    element.setNodeId(3,n3);
    model.addElement(element,true,groupID);
}

//! Add triangle to surface group.
void addTriangle(RModel &model, uint groupID, uint n0, uint n1, uint n2)
{
    RElement element(R_ELEMENT_TRI1);
    element.setNodeId(0,n0);
    element.setNodeId(1,n1);
    element.setNodeId(2,n2);
    model.addElement(element,true,groupID);
}

//! Build a unit cube of tetrahedra (each cell split into six) filled with gas and
//! a uniform particle rate. No velocity field is present, so the concentration
//! stays uniform and each node obeys dc/dt = f * (1 - c/cSat).
//! The box has length 1 with nx cells along X and width with nyz cells along Y and Z.
//! Tetrahedra only - hexahedra and wedges have no shape functions yet.
RModel buildCube(double particleRate, double maximumSaturation, bool sourceInFirstSlice = false, uint nx = nDivisions, uint nyz = nDivisions, double width = 1.0)
{
    auto nodeIndex = [nx,nyz](uint i, uint j, uint k) { return i + (nx+1)*(j + (nyz+1)*k); };

    RModel model;

    model.setNNodes((nx+1)*(nyz+1)*(nyz+1));
    for (uint k=0;k<=nyz;k++)
    {
        for (uint j=0;j<=nyz;j++)
        {
            for (uint i=0;i<=nx;i++)
            {
                model.getNode(nodeIndex(i,j,k)).set(double(i)/double(nx),
                                                    width*double(j)/double(nyz),
                                                    width*double(k)/double(nyz));
            }
        }
    }

    for (uint k=0;k<nyz;k++)
    {
        for (uint j=0;j<nyz;j++)
        {
            for (uint i=0;i<nx;i++)
            {
                uint h[8] = { nodeIndex(i  ,j  ,k  ),
                              nodeIndex(i+1,j  ,k  ),
                              nodeIndex(i+1,j+1,k  ),
                              nodeIndex(i  ,j+1,k  ),
                              nodeIndex(i  ,j  ,k+1),
                              nodeIndex(i+1,j  ,k+1),
                              nodeIndex(i+1,j+1,k+1),
                              nodeIndex(i  ,j+1,k+1) };

                // Cells beyond the first slice in X form a second volume without the source.
                uint groupID = (sourceInFirstSlice && i > 0) ? 1 : 0;

                // Split around the main diagonal 0-6.
                addTetra(model,groupID,h[0],h[1],h[2],h[6]);
                addTetra(model,groupID,h[0],h[2],h[3],h[6]);
                addTetra(model,groupID,h[0],h[3],h[7],h[6]);
                addTetra(model,groupID,h[0],h[7],h[4],h[6]);
                addTetra(model,groupID,h[0],h[4],h[5],h[6]);
                addTetra(model,groupID,h[0],h[5],h[1],h[6]);
            }
        }
    }

    RMaterial material(RMaterial::Gas);
    material.setName("Air");

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);

    RVolume &volume = model.getVolume(0);
    volume.setName("Room");
    volume.setMaterial(material);

    RBoundaryCondition rateBc(R_BOUNDARY_CONDITION_PARTICLE_RATE);
    rateBc.getComponent(rateBc.findComponentPosition(R_VARIABLE_PARTICLE_RATE)).add(0.0,particleRate);
    volume.addBoundaryCondition(rateBc);

    if (sourceInFirstSlice)
    {
        RVolume &cleanVolume = model.getVolume(1);
        cleanVolume.setName("Clean");
        cleanVolume.setMaterial(material);
    }

    model.getProblemSetup().getParticleSetup().setMaximumSaturation(maximumSaturation);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_FLUID_PARTICLE));

    model.getMatrixSolverConf(RMatrixSolverConf::GMRES).setOutputFrequency(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    return model;
}

//! Add triangles on the faces X = 0 and X = 1 of the cube from buildCube() as surfaces
//! with prescribed particle concentration, or free (natural) when prescribeConcentration is false.
void addEndFaces(RModel &model, uint nx, uint nyz, double cStart, double cEnd, bool prescribeConcentration = true)
{
    auto nodeIndex = [nx,nyz](uint i, uint j, uint k) { return i + (nx+1)*(j + (nyz+1)*k); };

    for (uint face=0;face<2;face++)
    {
        uint i = (face == 0) ? 0 : nx;
        for (uint k=0;k<nyz;k++)
        {
            for (uint j=0;j<nyz;j++)
            {
                addTriangle(model,face,nodeIndex(i,j,k),nodeIndex(i,j+1,k),nodeIndex(i,j+1,k+1));
                addTriangle(model,face,nodeIndex(i,j,k),nodeIndex(i,j+1,k+1),nodeIndex(i,j,k+1));
            }
        }

        RSurface &surface = model.getSurface(face);
        surface.setName(face == 0 ? "Start" : "End");

        if (prescribeConcentration)
        {
            RBoundaryCondition concentrationBc(R_BOUNDARY_CONDITION_PARTICLE_CONCENTRATION);
            concentrationBc.getComponent(concentrationBc.findComponentPosition(R_VARIABLE_PARTICLE_CONCENTRATION)).add(0.0,face == 0 ? cStart : cEnd);
            surface.addBoundaryCondition(concentrationBc);
        }
    }
}

void setTransient(RModel &model, uint nTimeSteps, double timeStepSize)
{
    RTimeSolver &timeSolver = model.getTimeSolver();
    timeSolver.setEnabled(true);
    timeSolver.setTimeMarchApproximation(R_TIME_MARCH_CENTRAL);
    timeSolver.setTimes(RTimeSolver::findTimesVector(nTimeSteps,0.0,timeStepSize));
    timeSolver.setInputNTimeSteps(nTimeSteps);
    timeSolver.setInputStartTime(0.0);
    timeSolver.setInputTimeStepSize(timeStepSize);
    timeSolver.setOutputFrequency(0);
}

//! Store node variable with given values.
void setNodeVariable(RModel &model, RVariableType variableType, const std::vector<RRVector> &values)
{
    uint position = model.addVariable(variableType);
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(uint(values.size()),model.getNNodes());
    for (uint i=0;i<values.size();i++)
    {
        for (uint j=0;j<model.getNNodes();j++)
        {
            variable.setValue(i,j,values[i][j]);
        }
    }
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

class TestSolverFluidParticle : public QObject
{
    Q_OBJECT

private slots:

    void steadySourceReachesSaturation();
    void transientSourceApproachesSaturation();
    void transientUnlimitedSourceIsLinear();
    void transientSmallCourantNumberIsStable();
    void transientNonSolenoidalFlowIsBounded();
    void transientBackflowIsBounded();
    void transientConvergingFlowCreatesNoMaximum();
    void transientBackflowAtNoSlipEdgeIsBounded();
    void steadyDiffusionIsLinear();
    void steadyFluxDiffusionIsLinear();
    void steadyMassTransferDiffusionIsLinear();
    void steadyAdvectionDiffusionProfile_data();
    void steadyAdvectionDiffusionProfile();
};


//! Without transport the steady balance f * (1 - c/cSat) = 0 holds only at saturation.
void TestSolverFluidParticle::steadySourceReachesSaturation()
{
    const double particleRate = 1.0e-3; // [kg/(m^3*s)]
    const double maximumSaturation = 0.05; // [kg/m^3]

    RModel model = buildCube(particleRate,maximumSaturation);
    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    QVERIFY(model.findVariable(R_VARIABLE_RELATIVE_SATURATION) != RConstants::eod);

    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double s = nodeValue(model,R_VARIABLE_RELATIVE_SATURATION,i);
        QVERIFY2(std::fabs(c-maximumSaturation) < 1.0e-6*maximumSaturation,
                 qPrintable(QString("Node %1: c = %2, expected %3").arg(i).arg(c).arg(maximumSaturation)));
        QVERIFY2(std::fabs(s-1.0) < 1.0e-6,
                 qPrintable(QString("Node %1: relative saturation = %2, expected 1").arg(i).arg(s)));
    }
}


//! Uniform concentration follows the Crank-Nicolson recursion of dc/dt = f * (1 - c/cSat)
//! and never exceeds saturation.
void TestSolverFluidParticle::transientSourceApproachesSaturation()
{
    const double particleRate = 1.0e-2; // [kg/(m^3*s)]
    const double maximumSaturation = 0.05; // [kg/m^3]
    const uint nTimeSteps = 20;
    const double dt = 1.0; // [s]

    RModel model = buildCube(particleRate,maximumSaturation);
    setTransient(model,nTimeSteps,dt);

    RSolver solver(model,QString(),QString());
    solver.run();

    double lambda = dt * particleRate / maximumSaturation;
    double expected = 0.0;
    for (uint n=0;n<nTimeSteps;n++)
    {
        expected = (dt * particleRate + expected * (1.0 - 0.5*lambda)) / (1.0 + 0.5*lambda);
    }
    // Sanity - after f*t/cSat = 4 the continuous solution is at 98 % of saturation.
    QVERIFY(expected > 0.95*maximumSaturation && expected < maximumSaturation);

    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double s = nodeValue(model,R_VARIABLE_RELATIVE_SATURATION,i);
        QVERIFY2(std::fabs(c-expected) < 1.0e-6*maximumSaturation,
                 qPrintable(QString("Node %1: c = %2, expected %3").arg(i).arg(c).arg(expected)));
        QVERIFY2(c <= maximumSaturation,
                 qPrintable(QString("Node %1: c = %2 exceeds saturation %3").arg(i).arg(c).arg(maximumSaturation)));
        QVERIFY2(std::fabs(s-c/maximumSaturation) < 1.0e-12,
                 qPrintable(QString("Node %1: relative saturation = %2, expected %3").arg(i).arg(s).arg(c/maximumSaturation)));
    }
}

//! Zero maximum saturation keeps the former unlimited behaviour, c = f * t.
void TestSolverFluidParticle::transientUnlimitedSourceIsLinear()
{
    const double particleRate = 1.0e-2; // [kg/(m^3*s)]
    const uint nTimeSteps = 20;
    const double dt = 1.0; // [s]

    RModel model = buildCube(particleRate,0.0);
    setTransient(model,nTimeSteps,dt);

    RSolver solver(model,QString(),QString());
    solver.run();

    QVERIFY(model.findVariable(R_VARIABLE_RELATIVE_SATURATION) == RConstants::eod);

    double expected = particleRate * dt * double(nTimeSteps);
    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        QVERIFY2(std::fabs(c-expected) < 1.0e-6*expected,
                 qPrintable(QString("Node %1: c = %2, expected %3").arg(i).arg(c).arg(expected)));
    }
}

//! In a nearly stagnant, irregular flow (element Courant number ~1e-5) transport is negligible,
//! so the concentration must match the solution without any flow. SUPG parameter which is not
//! limited by the time step (tau = h/(2|v|)) keeps the SUPG weighted mass term at full strength
//! with a direction changing from node to node, which spreads the contaminant spuriously.
void TestSolverFluidParticle::transientSmallCourantNumberIsStable()
{
    const double particleRate = 1.0e-2; // [kg/(m^3*s)]
    const uint nTimeSteps = 50;
    const double dt = 0.01; // [s]
    const double speed = 1.0e-3; // [m/s]

    std::vector<double> cReference;
    std::vector<double> cFlow;

    for (double s : { 0.0, speed })
    {
        RModel model = buildCube(particleRate,0.0,true);
        setTransient(model,nTimeSteps,dt);

        // Solver clears results on construction, velocity has to be set afterwards.
        RSolver solver(model,QString(),QString());

        // Velocity of constant magnitude with direction changing from node to node, as in a recirculation zone.
        std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
        for (uint i=0;i<model.getNNodes();i++)
        {
            RR3Vector direction(std::sin(1.7*double(i)+0.3),std::cos(2.3*double(i)+1.1),std::sin(3.1*double(i)+2.0));
            direction.normalize();
            for (uint j=0;j<3;j++)
            {
                velocity[j][i] = s * direction[j];
            }
        }
        setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

        solver.run();

        std::vector<double> &c = (s == 0.0) ? cReference : cFlow;
        for (uint i=0;i<model.getNNodes();i++)
        {
            c.push_back(nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i));
        }
    }

    double cMax = particleRate * dt * double(nTimeSteps);
    // Sanity - the source fills the first slice while the rest stays clean.
    QVERIFY(*std::max_element(cReference.begin(),cReference.end()) > 0.5*cMax);
    QVERIFY(*std::min_element(cReference.begin(),cReference.end()) < 0.1*cMax);

    double maxDifference = 0.0;
    for (uint i=0;i<cFlow.size();i++)
    {
        QVERIFY2(std::isfinite(cFlow[i]),qPrintable(QString("Node %1: non finite concentration").arg(i)));
        maxDifference = std::max(maxDifference,std::fabs(cFlow[i]-cReference[i]));
    }
    QVERIFY2(maxDifference < 0.01*cMax,
             qPrintable(QString("Concentration differs from no flow solution by %1 (%2 %)").arg(maxDifference).arg(100.0*maxDifference/cMax)));
}

//! Pure diffusion between two faces with prescribed concentration gives linear
//! profile, which linear elements reproduce exactly.
void TestSolverFluidParticle::steadyDiffusionIsLinear()
{
    const uint nx = 6;
    const uint nyz = 2;

    RModel model = buildCube(0.0,0.0,false,nx,nyz);
    addEndFaces(model,nx,nyz,0.0,1.0);
    model.getProblemSetup().getParticleSetup().setDiffusionCoefficient(1.0e-3);
    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double x = model.getNode(i).getX();
        QVERIFY2(std::fabs(c-x) < 1.0e-6,
                 qPrintable(QString("Node %1 (x = %2): c = %3, expected %4").arg(i).arg(x).arg(c).arg(x)));
    }
}

//! Pure diffusion from a face with prescribed particle flux q to a face with zero
//! concentration gives linear profile c = q/D * (1 - x).
void TestSolverFluidParticle::steadyFluxDiffusionIsLinear()
{
    const uint nx = 6;
    const uint nyz = 2;
    const double diffusionCoefficient = 1.0e-3;
    const double flux = 2.0e-3;

    RModel model = buildCube(0.0,0.0,false,nx,nyz);
    addEndFaces(model,nx,nyz,0.0,0.0,false);

    RBoundaryCondition fluxBc(R_BOUNDARY_CONDITION_PARTICLE_FLUX);
    fluxBc.getComponent(fluxBc.findComponentPosition(R_VARIABLE_PARTICLE_FLUX)).add(0.0,flux);
    model.getSurface(0).addBoundaryCondition(fluxBc);

    RBoundaryCondition concentrationBc(R_BOUNDARY_CONDITION_PARTICLE_CONCENTRATION);
    concentrationBc.getComponent(concentrationBc.findComponentPosition(R_VARIABLE_PARTICLE_CONCENTRATION)).add(0.0,0.0);
    model.getSurface(1).addBoundaryCondition(concentrationBc);

    model.getProblemSetup().getParticleSetup().setDiffusionCoefficient(diffusionCoefficient);
    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double x = model.getNode(i).getX();
        double expected = flux / diffusionCoefficient * (1.0 - x);
        QVERIFY2(std::fabs(c-expected) < 1.0e-6,
                 qPrintable(QString("Node %1 (x = %2): c = %3, expected %4").arg(i).arg(x).arg(c).arg(expected)));
    }
}

//! Pure diffusion from a face with mass transfer beta * (cInf - c) to a face with zero
//! concentration gives linear profile c = c0 * (1 - x), where D * c0 = beta * (cInf - c0).
void TestSolverFluidParticle::steadyMassTransferDiffusionIsLinear()
{
    const uint nx = 6;
    const uint nyz = 2;
    const double diffusionCoefficient = 1.0e-3;
    const double massTransferCoefficient = 3.0e-3;
    const double ambientConcentration = 2.0;

    RModel model = buildCube(0.0,0.0,false,nx,nyz);
    addEndFaces(model,nx,nyz,0.0,0.0,false);

    RBoundaryCondition massTransferBc(R_BOUNDARY_CONDITION_PARTICLE_MASS_TRANSFER);
    massTransferBc.getComponent(massTransferBc.findComponentPosition(R_VARIABLE_PARTICLE_MASS_TRANSFER_COEFFICIENT)).add(0.0,massTransferCoefficient);
    massTransferBc.getComponent(massTransferBc.findComponentPosition(R_VARIABLE_PARTICLE_AMBIENT_CONCENTRATION)).add(0.0,ambientConcentration);
    model.getSurface(0).addBoundaryCondition(massTransferBc);

    RBoundaryCondition concentrationBc(R_BOUNDARY_CONDITION_PARTICLE_CONCENTRATION);
    concentrationBc.getComponent(concentrationBc.findComponentPosition(R_VARIABLE_PARTICLE_CONCENTRATION)).add(0.0,0.0);
    model.getSurface(1).addBoundaryCondition(concentrationBc);

    model.getProblemSetup().getParticleSetup().setDiffusionCoefficient(diffusionCoefficient);
    model.getTimeSolver().setEnabled(false);

    RSolver solver(model,QString(),QString());
    solver.run();

    double c0 = massTransferCoefficient * ambientConcentration / (diffusionCoefficient + massTransferCoefficient);

    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double x = model.getNode(i).getX();
        double expected = c0 * (1.0 - x);
        QVERIFY2(std::fabs(c-expected) < 1.0e-6,
                 qPrintable(QString("Node %1 (x = %2): c = %3, expected %4").arg(i).arg(x).arg(c).arg(expected)));
    }
}

void TestSolverFluidParticle::steadyAdvectionDiffusionProfile_data()
{
    QTest::addColumn<double>("pecletNumber");
    QTest::addColumn<double>("tolerance");
    // Diffusion dominated - SUPG is nearly off.
    QTest::newRow("Pe = 1") << 1.0 << 1.0e-5;
    // Advection dominated with boundary layer at the end - SUPG and diffusion combined.
    QTest::newRow("Pe = 20") << 20.0 << 1.0e-2;
}

//! Uniform flow along X with diffusion between faces c(0) = 0 and c(1) = 1 has
//! profile c = (exp(Pe*x) - 1) / (exp(Pe) - 1), Pe = |v| L / D.
void TestSolverFluidParticle::steadyAdvectionDiffusionProfile()
{
    QFETCH(double,pecletNumber);
    QFETCH(double,tolerance);

    const uint nx = 40;
    const uint nyz = 1;
    const double speed = 1.0; // [m/s]

    // Cross section of one cubic cell - flat tetrahedra would couple the corners and spoil the 1D solution.
    RModel model = buildCube(0.0,0.0,false,nx,nyz,1.0/double(nx));
    addEndFaces(model,nx,nyz,0.0,1.0);
    model.getProblemSetup().getParticleSetup().setDiffusionCoefficient(speed / pecletNumber);
    // Steady diffusion converges slowly - default GMRES limit (10 x 10 iterations) is not enough.
    model.getMatrixSolverConf(RMatrixSolverConf::GMRES).setNOuterIterations(1000);
    model.getTimeSolver().setEnabled(false);

    // Solver clears results on construction, velocity has to be set afterwards.
    RSolver solver(model,QString(),QString());

    std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
    velocity[0].fill(speed);
    setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

    solver.run();

    double maxError = 0.0;
    for (uint i=0;i<model.getNNodes();i++)
    {
        double c = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        double x = model.getNode(i).getX();
        double expected = std::expm1(pecletNumber*x) / std::expm1(pecletNumber);
        QVERIFY2(std::isfinite(c),qPrintable(QString("Node %1: non finite concentration").arg(i)));
        maxError = std::max(maxError,std::fabs(c-expected));
    }
    QVERIFY2(maxError < tolerance,
             qPrintable(QString("Maximum error %1 exceeds %2").arg(maxError).arg(tolerance)));
}

//! Velocity of a stabilised flow solution is not exactly divergence free. Pure advection on it
//! must not grow on its own: with no flow through the walls the solution may overshoot the
//! bound f * t locally (the scheme is not monotone), but doubling the run must not do more than
//! about double the peak. Advection term which is not energy neutral for such a velocity grows
//! exponentially once the SUPG diffusion is small (small time step).
void TestSolverFluidParticle::transientNonSolenoidalFlowIsBounded()
{
    const double particleRate = 1.0; // [kg/(m^3*s)]
    const uint nTimeSteps = 200;
    const double dt = 0.01; // [s]
    const double speed = 1.0; // [m/s]
    const uint n = 6;

    std::vector<double> cMax;

    for (uint nSteps : { nTimeSteps, 2*nTimeSteps })
    {
        RModel model = buildCube(particleRate,0.0,true,n,n);
        setTransient(model,nSteps,dt);

        // Solver clears results on construction, velocity has to be set afterwards.
        RSolver solver(model,QString(),QString());

        // Direction changing from node to node - far from divergence free - and zero on the walls
        // so that nothing flows in or out.
        std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
        for (uint i=0;i<model.getNNodes();i++)
        {
            const RNode &node = model.getNode(i);
            bool onWall = false;
            for (double coordinate : { node.getX(), node.getY(), node.getZ() })
            {
                onWall = onWall || coordinate < 1.0e-9 || coordinate > 1.0 - 1.0e-9;
            }
            if (onWall)
            {
                continue;
            }
            RR3Vector direction(std::sin(1.7*double(i)+0.3),std::cos(2.3*double(i)+1.1),std::sin(3.1*double(i)+2.0));
            direction.normalize();
            for (uint j=0;j<3;j++)
            {
                velocity[j][i] = speed * direction[j];
            }
        }
        setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

        solver.run();

        double c = 0.0;
        for (uint i=0;i<model.getNNodes();i++)
        {
            double ci = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
            QVERIFY2(std::isfinite(ci),qPrintable(QString("Node %1: non finite concentration").arg(i)));
            c = std::max(c,std::fabs(ci));
        }
        cMax.push_back(c);
    }

    double cBound = particleRate * dt * double(2*nTimeSteps);
    QVERIFY2(cMax[1] < 3.0*cMax[0],
             qPrintable(QString("Peak concentration grows faster than the injected mass: %1 after %2 steps, %3 after %4 steps")
                        .arg(cMax[0]).arg(nTimeSteps).arg(cMax[1]).arg(2*nTimeSteps)));
    QVERIFY2(cMax[1] < 5.0*cBound,
             qPrintable(QString("Peak concentration %1 far beyond the bound f * t = %2").arg(cMax[1]).arg(cBound)));
}

//! Flow entering through a free boundary (backflow at a pressure outlet) carries no prescribed
//! concentration. Without an inflow treatment the advection boundary term feeds energy into the
//! solution there and it grows without bound. Here all flow enters through the free face X = 1
//! and leaves through the free face X = 0; the source sits next to the exit, so the steady
//! concentration is about f * h / |v| and the field must stay near it.
void TestSolverFluidParticle::transientBackflowIsBounded()
{
    const double particleRate = 1.0; // [kg/(m^3*s)]
    const uint nTimeSteps = 200;
    const double dt = 0.01; // [s]
    const double speed = 1.0; // [m/s]
    const uint n = 6;

    std::vector<double> cMax;

    for (uint nSteps : { nTimeSteps, 2*nTimeSteps })
    {
        RModel model = buildCube(particleRate,0.0,true,n,n);
        addEndFaces(model,n,n,0.0,0.0,false);
        setTransient(model,nSteps,dt);

        // Solver clears results on construction, velocity has to be set afterwards.
        RSolver solver(model,QString(),QString());

        std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
        velocity[0].fill(-speed);
        setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

        solver.run();

        double c = 0.0;
        for (uint i=0;i<model.getNNodes();i++)
        {
            double ci = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
            QVERIFY2(std::isfinite(ci),qPrintable(QString("Node %1: non finite concentration").arg(i)));
            c = std::max(c,std::fabs(ci));
        }
        cMax.push_back(c);
    }

    // Residence time in the source slice is (1/n) / |v|.
    double cSteady = particleRate / (double(n) * speed);
    QVERIFY2(cMax[1] < 3.0*cSteady && cMax[1] < 1.5*cMax[0],
             qPrintable(QString("Concentration grows: %1 after %2 steps, %3 after %4 steps, steady value about %5")
                        .arg(cMax[0]).arg(nTimeSteps).arg(cMax[1]).arg(2*nTimeSteps).arg(cSteady)));
}

//! In converging flow (div(v) < 0) a uniform concentration is the exact solution of the advection
//! equation. The skew-symmetric advection term 0.5 * div(v) * C acts there as a source and the
//! concentration grows, piling up where a discrete flow field is far from divergence free. Here
//! the flow converges to the cube centre and is zero on the walls. Next to the walls the flow
//! expands and the divergence sink lowers the concentration; the gradient it creates is then
//! steepened by the converging flow and Galerkin advection overshoots - a bounded overshoot (about
//! 2) is expected, the source of the skew-symmetric form drives it to about 7.
void TestSolverFluidParticle::transientConvergingFlowCreatesNoMaximum()
{
    const uint nTimeSteps = 100;
    const double dt = 0.01; // [s]
    const double rate = 5.0; // [1/s] convergence rate, div(v) = -3 * rate
    const uint n = 6;

    RModel model = buildCube(0.0,0.0,false,n,n);

    RInitialCondition concentrationIc(R_INITIAL_CONDITION_PARTICLE_CONCENTRATION);
    concentrationIc.getComponent(concentrationIc.findComponentPosition(R_VARIABLE_PARTICLE_CONCENTRATION)).add(0.0,1.0);
    model.getVolume(0).addInitialCondition(concentrationIc);

    setTransient(model,nTimeSteps,dt);

    // Solver clears results on construction, velocity has to be set afterwards.
    RSolver solver(model,QString(),QString());

    std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
    for (uint i=0;i<model.getNNodes();i++)
    {
        const RNode &node = model.getNode(i);
        bool onWall = false;
        for (double coordinate : { node.getX(), node.getY(), node.getZ() })
        {
            onWall = onWall || coordinate < 1.0e-9 || coordinate > 1.0 - 1.0e-9;
        }
        if (onWall)
        {
            continue;
        }
        velocity[0][i] = -rate * (node.getX() - 0.5);
        velocity[1][i] = -rate * (node.getY() - 0.5);
        velocity[2][i] = -rate * (node.getZ() - 0.5);
    }
    setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

    solver.run();

    double cMax = 0.0;
    for (uint i=0;i<model.getNNodes();i++)
    {
        double ci = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
        QVERIFY2(std::isfinite(ci),qPrintable(QString("Node %1: non finite concentration").arg(i)));
        cMax = std::max(cMax,ci);
    }
    QVERIFY2(cMax < 3.0,
             qPrintable(QString("Converging flow creates a maximum: %1, initial value 1").arg(cMax)));
}

//! Flow enters through the free face X = 1 with a profile that is zero on the face edges, as at a
//! no-slip wall bordering an outlet with backflow. The energy fed in by the advection boundary
//! term over a face reaches its edge nodes although their own velocity is zero, so the inflow
//! treatment must damp them as well. Clean fluid entering through X = 1 flushes the initial
//! concentration out through X = 0; no value may grow beyond the initial one.
void TestSolverFluidParticle::transientBackflowAtNoSlipEdgeIsBounded()
{
    const uint nTimeSteps = 200;
    const double dt = 0.01; // [s]
    const double speed = 2.0; // [m/s]
    const uint n = 6;

    std::vector<double> cMax;

    for (uint nSteps : { nTimeSteps, 2*nTimeSteps })
    {
        RModel model = buildCube(0.0,0.0,false,n,n);
        addEndFaces(model,n,n,0.0,0.0,false);

        RInitialCondition concentrationIc(R_INITIAL_CONDITION_PARTICLE_CONCENTRATION);
        concentrationIc.getComponent(concentrationIc.findComponentPosition(R_VARIABLE_PARTICLE_CONCENTRATION)).add(0.0,1.0);
        model.getVolume(0).addInitialCondition(concentrationIc);

        setTransient(model,nSteps,dt);

        // Solver clears results on construction, velocity has to be set afterwards.
        RSolver solver(model,QString(),QString());

        // Divergence free profile, entering through X = 1 and leaving through X = 0, zero on walls Y and Z.
        std::vector<RRVector> velocity(3,RRVector(model.getNNodes(),0.0));
        for (uint i=0;i<model.getNNodes();i++)
        {
            const RNode &node = model.getNode(i);
            velocity[0][i] = -speed * std::sin(M_PI*node.getY()) * std::sin(M_PI*node.getZ());
        }
        setNodeVariable(model,R_VARIABLE_VELOCITY,velocity);

        solver.run();

        double c = 0.0;
        for (uint i=0;i<model.getNNodes();i++)
        {
            double ci = nodeValue(model,R_VARIABLE_PARTICLE_CONCENTRATION,i);
            QVERIFY2(std::isfinite(ci),qPrintable(QString("Node %1: non finite concentration").arg(i)));
            c = std::max(c,std::fabs(ci));
        }
        cMax.push_back(c);
    }

    qInfo("Peak concentration %g after %u steps, %g after %u steps",cMax[0],nTimeSteps,cMax[1],2*nTimeSteps);

    // Galerkin advection is not monotone - a bounded overshoot of the initial value is expected,
    // growth is not.
    QVERIFY2(cMax[1] < 1.2*cMax[0] && cMax[1] < 3.0,
             qPrintable(QString("Concentration grows at the inflow edge: %1 after %2 steps, %3 after %4 steps, initial value 1")
                        .arg(cMax[0]).arg(nTimeSteps).arg(cMax[1]).arg(2*nTimeSteps)));
}

QTEST_APPLESS_MAIN(TestSolverFluidParticle)

#include "tst_solver_fluid_particle.moc"

#include <QtTest>

#include <cmath>

#include "rml_model.h"
#include "rsolver.h"

namespace
{

// A solid slab and a fluid slab side by side along the x axis, meshed
// conformally, with the Forced convection condition on the interface between
// them. The heat flow is one dimensional, which gives the coupled solution in
// closed form.
constexpr double solidLength = 0.1;          // [m]
constexpr double fluidLength = 0.1;          // [m]
constexpr double side = 0.02;                // [m]
constexpr uint   nSolidSlices = 10;
constexpr uint   nFluidSlices = 40;
constexpr uint   nSlices = nSolidSlices + nFluidSlices;

constexpr double solidConductivity = 2.0;    // [W/(m*K)]
constexpr double fluidConductivity = 0.6;    // [W/(m*K)]
constexpr double fluidDensity = 1000.0;      // [kg/m^3]
constexpr double fluidCapacity = 4180.0;     // [J/(kg*K)]

constexpr double hotTemperature = 400.0;     // [K] at x = 0, the solid end
constexpr double coldTemperature = 300.0;    // [K] at x = L, the fluid end

uint nodeId(uint i, uint j, uint k)
{
    return 4*i + 2*j + k;
}

double sliceX(uint i)
{
    if (i <= nSolidSlices)
    {
        return solidLength * double(i) / double(nSolidSlices);
    }
    return solidLength + fluidLength * double(i - nSolidSlices) / double(nFluidSlices);
}

void addTetrahedron(RModel &model, uint groupID, uint n1, uint n2, uint n3, uint n4)
{
    RR3Vector a(model.getNode(n1).toVector());
    RR3Vector b(model.getNode(n2).toVector());
    RR3Vector c(model.getNode(n3).toVector());
    RR3Vector d(model.getNode(n4).toVector());

    RR3Vector ab,ac,ad,cross;
    RR3Vector::subtract(b,a,ab);
    RR3Vector::subtract(c,a,ac);
    RR3Vector::subtract(d,a,ad);
    RR3Vector::cross(ab,ac,cross);

    if (RR3Vector::dot(cross,ad) < 0.0)
    {
        std::swap(n3,n4);
    }

    RElement element(R_ELEMENT_TETRA1);
    element.setNodeId(0,n1);
    element.setNodeId(1,n2);
    element.setNodeId(2,n3);
    element.setNodeId(3,n4);
    model.addElement(element,true,groupID);
}

//! Add the two triangles of the square face at slice i. The diagonal is the
//! one the Kuhn decomposition of the neighbouring bricks cuts the face along.
void addFace(RModel &model, uint groupID, uint i)
{
    RElement t1(R_ELEMENT_TRI1);
    t1.setNodeId(0,nodeId(i,0,0));
    t1.setNodeId(1,nodeId(i,1,0));
    t1.setNodeId(2,nodeId(i,1,1));
    model.addElement(t1,true,groupID);

    RElement t2(R_ELEMENT_TRI1);
    t2.setNodeId(0,nodeId(i,0,0));
    t2.setNodeId(1,nodeId(i,0,1));
    t2.setNodeId(2,nodeId(i,1,1));
    model.addElement(t2,true,groupID);
}

RMaterialProperty property(RMaterialProperty::Type type, double value)
{
    RMaterialProperty materialProperty(type);
    materialProperty.clear();
    materialProperty.add(293.15,value);
    return materialProperty;
}

RBoundaryCondition temperatureBc(double temperature)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_TEMPERATURE);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_TEMPERATURE)).add(0.0,temperature);
    return bc;
}

//! Build the slabs. Volume 0 is the solid, volume 1 the fluid. Surface 0 is
//! the interface, surface 1 the hot solid end and surface 2 the cold fluid end.
RModel buildSlabs(bool forcedConvection, bool fluidWithEmissivity)
{
    RModel model;

    model.setNNodes(4*(nSlices+1));
    for (uint i=0;i<=nSlices;i++)
    {
        for (uint j=0;j<2;j++)
        {
            for (uint k=0;k<2;k++)
            {
                model.getNode(nodeId(i,j,k)).set(sliceX(i),double(j)*side,double(k)*side);
            }
        }
    }

    for (uint i=0;i<nSlices;i++)
    {
        uint groupID = (i < nSolidSlices) ? 0 : 1;

        uint v000 = nodeId(i,0,0);
        uint v001 = nodeId(i,0,1);
        uint v010 = nodeId(i,1,0);
        uint v011 = nodeId(i,1,1);
        uint v100 = nodeId(i+1,0,0);
        uint v101 = nodeId(i+1,0,1);
        uint v110 = nodeId(i+1,1,0);
        uint v111 = nodeId(i+1,1,1);

        addTetrahedron(model,groupID,v000,v100,v110,v111);
        addTetrahedron(model,groupID,v000,v110,v010,v111);
        addTetrahedron(model,groupID,v000,v010,v011,v111);
        addTetrahedron(model,groupID,v000,v011,v001,v111);
        addTetrahedron(model,groupID,v000,v001,v101,v111);
        addTetrahedron(model,groupID,v000,v101,v100,v111);
    }

    addFace(model,0,nSolidSlices);
    addFace(model,1,0);
    addFace(model,2,nSlices);

    RMaterial solid(RMaterial::Solid);
    solid.setName("Solid");
    solid.add(property(RMaterialProperty::Density,2700.0));
    solid.add(property(RMaterialProperty::HeatCapacity,900.0));
    solid.add(property(RMaterialProperty::ThermalConductivity,solidConductivity));
    solid.add(property(RMaterialProperty::Emissivity,0.5));
    model.getVolume(0).setName("Solid");
    model.getVolume(0).setMaterial(solid);

    // The state is left unspecified, as in the shipped material database - the
    // dynamic viscosity alone has to make it a fluid.
    RMaterial fluid;
    fluid.setName("Fluid");
    fluid.add(property(RMaterialProperty::Density,fluidDensity));
    fluid.add(property(RMaterialProperty::DynamicViscosity,1.0e-3));
    fluid.add(property(RMaterialProperty::HeatCapacity,fluidCapacity));
    fluid.add(property(RMaterialProperty::ThermalConductivity,fluidConductivity));
    if (fluidWithEmissivity)
    {
        // Carries every property the heat solver asks for, as mercury does.
        fluid.add(property(RMaterialProperty::Emissivity,0.1));
    }
    model.getVolume(1).setName("Fluid");
    model.getVolume(1).setMaterial(fluid);

    model.getSurface(0).setName("Interface");
    if (forcedConvection)
    {
        model.getSurface(0).addBoundaryCondition(RBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED));
    }

    model.getSurface(1).setName("Hot end");
    model.getSurface(1).addBoundaryCondition(temperatureBc(hotTemperature));

    model.getSurface(2).setName("Cold end");
    model.getSurface(2).addBoundaryCondition(temperatureBc(coldTemperature));

    for (RMatrixSolverConf::Type type : {RMatrixSolverConf::CG, RMatrixSolverConf::GMRES})
    {
        model.getMatrixSolverConf(type).setNInnerIterations(500);
        model.getMatrixSolverConf(type).setNOuterIterations(5000);
        model.getMatrixSolverConf(type).setSolverCvgValue(1.0e-14);
        model.getMatrixSolverConf(type).setOutputFrequency(0);
    }

    model.getTimeSolver().setEnabled(false);

    return model;
}

//! Prescribe a uniform flow along x in the fluid, as the flow task would have
//! stored it. The run is marked as a restart, otherwise the solver clears all
//! results - the prescribed velocity with them.
void setVelocity(RModel &model, double velocity)
{
    uint position = model.addVariable(R_VARIABLE_VELOCITY);
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(3,model.getNNodes());
    for (uint i=0;i<model.getNNodes();i++)
    {
        variable.setValue(0,i,(i >= nodeId(nSolidSlices,0,0)) ? velocity : 0.0);
        variable.setValue(1,i,0.0);
        variable.setValue(2,i,0.0);
    }
    model.getProblemSetup().setRestart(true);
}

void setCoupledTasks(RModel &model, uint nIterations)
{
    RProblemTaskItem root;
    root.setNIterations(nIterations);
    root.setCvgValue(1.0e-9);
    root.addChild(RProblemTaskItem(R_PROBLEM_FLUID_HEAT));
    root.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    model.setProblemTaskTree(root);
}

double nodeTemperature(const RModel &model, uint nodeID)
{
    uint position = model.findVariable(R_VARIABLE_TEMPERATURE);
    return model.getVariable(position).getValue(0,nodeID);
}

//! Interface temperature of the coupled problem. The solid conducts, the fluid
//! advects and conducts - with a Peclet number of zero the fluid only conducts.
double exactInterfaceTemperature(double velocity)
{
    double solidConductance = solidConductivity / solidLength;
    double pe = fluidDensity * fluidCapacity * velocity * fluidLength / fluidConductivity;
    double fluidConductance = (std::fabs(pe) < 1.0e-12)
                            ? fluidConductivity / fluidLength
                            : fluidConductivity * pe / (fluidLength * (std::exp(pe) - 1.0));
    return (solidConductance * hotTemperature + fluidConductance * coldTemperature) / (solidConductance + fluidConductance);
}

bool runSolver(RModel &model, QString &errorMessage)
{
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        errorMessage = error.getMessage();
        return false;
    }
    return true;
}

} // namespace

class TestSolverHeatCoupling : public QObject
{
    Q_OBJECT

private slots:

    void heatSolvesSolidsOnly();
    void fluidHeatSource();
    void conductingFluid();
    void advectingFluid();
};

void TestSolverHeatCoupling::fluidHeatSource()
{
    // A fluid heat task alone, with the fluid at rest and heated uniformly. The
    // interface side is insulated and the far end held cold, so the temperature
    // is a parabola peaking at the interface. A heat source used to cool the
    // fluid, because the conduction entered the element matrix with the wrong
    // sign.
    constexpr double sourceDensity = 1200.0; // [W/m^3]

    RModel model = buildSlabs(false,false);
    setVelocity(model,0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_FLUID_HEAT));

    RBoundaryCondition heatBc(R_BOUNDARY_CONDITION_HEAT);
    heatBc.getComponent(heatBc.findComponentPosition(R_VARIABLE_HEAT)).add(0.0,sourceDensity*fluidLength*side*side);
    model.getVolume(1).addBoundaryCondition(heatBc);

    QString errorMessage;
    QVERIFY2(runSolver(model,errorMessage),errorMessage.toUtf8().constData());

    double expected = coldTemperature + sourceDensity * fluidLength * fluidLength / (2.0 * fluidConductivity);
    double computed = nodeTemperature(model,nodeId(nSolidSlices,0,0));
    QVERIFY2(std::fabs(computed - expected) < 0.01 * (expected - coldTemperature),
             qPrintable(QString("interface %1 K, expected %2 K").arg(computed).arg(expected)));
}

void TestSolverHeatCoupling::heatSolvesSolidsOnly()
{
    // A heat task alone, with no fluid heat task and an insulated interface.
    // The fluid carries every property the heat solver asks for, yet must not
    // conduct - the solid stays at the temperature of its hot end, and the fluid
    // is left where it started.
    RModel model = buildSlabs(false,true);
    setVelocity(model,0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_HEAT));

    QString errorMessage;
    QVERIFY2(runSolver(model,errorMessage),errorMessage.toUtf8().constData());

    for (uint i=0;i<=nSolidSlices;i++)
    {
        QVERIFY(std::fabs(nodeTemperature(model,nodeId(i,0,0)) - hotTemperature) < 1.0e-6);
    }
    double initTemperature = RVariable::getInitValue(R_VARIABLE_TEMPERATURE);
    for (uint i=nSolidSlices+1;i<nSlices;i++)
    {
        QCOMPARE(nodeTemperature(model,nodeId(i,1,1)),initTemperature);
    }
}

void TestSolverHeatCoupling::conductingFluid()
{
    // A fluid at rest is a plain conductor, and the coupled solution is the one
    // of two conductors in series - linear in each slab, which linear elements
    // reproduce exactly. The nodes of the interface see no velocity at all,
    // which left the former correlation with no heat transfer. Ten passes are
    // plenty with the relaxed wall temperature - plain alternation needs about
    // a hundred and sixty here, the solid conducting poorly against the first
    // fluid element.
    RModel model = buildSlabs(true,false);
    setVelocity(model,0.0);
    setCoupledTasks(model,10);

    QString errorMessage;
    QVERIFY2(runSolver(model,errorMessage),errorMessage.toUtf8().constData());

    double expected = exactInterfaceTemperature(0.0);
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            double computed = nodeTemperature(model,nodeId(nSolidSlices,j,k));
            QVERIFY2(std::fabs(computed - expected) < 1.0e-3,
                     qPrintable(QString("interface %1 K, expected %2 K").arg(computed).arg(expected)));
        }
    }

    // Linear in the solid.
    double midSolid = nodeTemperature(model,nodeId(nSolidSlices/2,0,0));
    QVERIFY(std::fabs(midSolid - 0.5*(hotTemperature + expected)) < 1.0e-3);
}

void TestSolverHeatCoupling::advectingFluid()
{
    // The cold fluid flows towards the interface and steepens the temperature
    // gradient in front of it, so the interface runs colder than with the fluid
    // at rest. The Peclet number of the fluid slab is 7, resolved by elements of
    // cell Peclet number well below one. The wall flux is the consistent one, so
    // the boundary layer costs no accuracy - a flux taken from the gradient of
    // the first fluid element misses the interface temperature by 1.5 K here.
    double velocity = -7.0 * fluidConductivity / (fluidDensity * fluidCapacity * fluidLength);

    RModel model = buildSlabs(true,false);
    setVelocity(model,velocity);
    setCoupledTasks(model,10);

    QString errorMessage;
    QVERIFY2(runSolver(model,errorMessage),errorMessage.toUtf8().constData());

    double expected = exactInterfaceTemperature(velocity);
    double atRest = exactInterfaceTemperature(0.0);
    QVERIFY(expected < atRest - 1.0);

    double computed = nodeTemperature(model,nodeId(nSolidSlices,0,0));
    QVERIFY2(std::fabs(computed - expected) < 1.0e-2,
             qPrintable(QString("interface %1 K, expected %2 K").arg(computed).arg(expected)));
}

QTEST_APPLESS_MAIN(TestSolverHeatCoupling)

#include "tst_solver_heat_coupling.moc"

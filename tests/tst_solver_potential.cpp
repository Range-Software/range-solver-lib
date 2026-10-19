#include <QtTest>

#include <cmath>

#include "rml_model.h"
#include "rsolver.h"

namespace
{

constexpr double columnLength = 2.0;    // [m]
constexpr double columnSide = 0.2;      // [m]
constexpr uint   nx = 10;               // divisions along the length
constexpr uint   nyz = 2;               // divisions across the side

constexpr double permeability = 1.0e-11; // [m^2]
constexpr double viscosity = 1.0e-3;     // [kg/(m*s)]
constexpr double density = 1000.0;       // [kg/m^3]
constexpr double conductivity = permeability / viscosity;

uint nodeIndex(uint i, uint j, uint k)
{
    return i + (nx+1)*(j + (nyz+1)*k);
}

//! Add tetrahedron with positive orientation.
void addTetra(RModel &model, uint n0, uint n1, uint n2, uint n3)
{
    RR3Vector p0(model.getNode(n0).toVector());
    RR3Vector p1(model.getNode(n1).toVector());
    RR3Vector p2(model.getNode(n2).toVector());
    RR3Vector p3(model.getNode(n3).toVector());

    RR3Vector a,b,c,cross;
    RR3Vector::subtract(p1,p0,a);
    RR3Vector::subtract(p2,p0,b);
    RR3Vector::subtract(p3,p0,c);
    RR3Vector::cross(a,b,cross);

    if (RR3Vector::dot(cross,c) < 0.0)
    {
        std::swap(n1,n2);
    }

    RElement element(R_ELEMENT_TETRA1);
    element.setNodeId(0,n0);
    element.setNodeId(1,n1);
    element.setNodeId(2,n2);
    element.setNodeId(3,n3);
    model.addElement(element,true,0);
}

void addTriangle(RModel &model, uint groupID, uint n0, uint n1, uint n2)
{
    RElement element(R_ELEMENT_TRI1);
    element.setNodeId(0,n0);
    element.setNodeId(1,n1);
    element.setNodeId(2,n2);
    model.addElement(element,true,groupID);
}

//! Build a porous column of tetrahedra along the X axis spanning
//! [0,L] x [0,a] x [0,a] with end faces X = 0 (surface 0) and X = L (surface 1).
RModel buildColumn()
{
    RModel model;

    model.setNNodes((nx+1)*(nyz+1)*(nyz+1));
    for (uint k=0;k<=nyz;k++)
    {
        for (uint j=0;j<=nyz;j++)
        {
            for (uint i=0;i<=nx;i++)
            {
                model.getNode(nodeIndex(i,j,k)).set(columnLength*double(i)/double(nx),
                                                    columnSide*double(j)/double(nyz),
                                                    columnSide*double(k)/double(nyz));
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

                addTetra(model,h[0],h[1],h[2],h[6]);
                addTetra(model,h[0],h[2],h[3],h[6]);
                addTetra(model,h[0],h[3],h[7],h[6]);
                addTetra(model,h[0],h[7],h[4],h[6]);
                addTetra(model,h[0],h[4],h[5],h[6]);
                addTetra(model,h[0],h[5],h[1],h[6]);
            }
        }
    }

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
    }
    model.getSurface(0).setName("Inlet");
    model.getSurface(1).setName("Outlet");

    RMaterial material(RMaterial::Solid);
    material.setName("Sand");

    RMaterialProperty permeabilityProperty(RMaterialProperty::PermeabilityToFluids);
    permeabilityProperty.add(293.15,permeability);
    material.add(permeabilityProperty);

    RMaterialProperty viscosityProperty(RMaterialProperty::DynamicViscosity);
    viscosityProperty.add(293.15,viscosity);
    material.add(viscosityProperty);

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);

    model.getVolume(0).setName("Column");
    model.getVolume(0).setMaterial(material);

    model.getTimeSolver().setEnabled(false);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_POTENTIAL));

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    return model;
}

//! Element groups carry the default gravity (0,0,-g) - replace it on the column.
void setGravity(RModel &model, double gx, double gy, double gz)
{
    REnvironmentCondition gravity(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_X)).add(0.0,gx);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y)).add(0.0,gy);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Z)).add(0.0,gz);
    model.getVolume(0).addEnvironmentCondition(gravity);
}

void setPressure(RModel &model, uint surfaceID, double pressure)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_PRESSURE)).add(0.0,pressure);
    model.getSurface(surfaceID).addBoundaryCondition(bc);
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
    catch (const std::exception &error)
    {
        errorMessage = QString(error.what());
        return false;
    }
    return true;
}

double nodePressure(const RModel &model, uint nodeID)
{
    uint position = model.findVariable(R_VARIABLE_PRESSURE);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(0,nodeID);
}

RR3Vector elementVelocity(const RModel &model, uint elementID)
{
    uint position = model.findVariable(R_VARIABLE_VELOCITY);
    if (position == RConstants::eod)
    {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return RR3Vector(nan,nan,nan);
    }
    const RVariable &variable = model.getVariable(position);
    return RR3Vector(variable.getValue(0,elementID),variable.getValue(1,elementID),variable.getValue(2,elementID));
}

//! Check that pressure is linear p(x) = p0 + (p1 - p0)*x/L and the velocity
//! of every volume element is v = (vx,0,0).
void verifyUniformFlow(const RModel &model, double p0, double p1, double vx)
{
    double pScale = std::max(std::abs(p0),std::abs(p1));
    for (uint i=0;i<model.getNNodes();i++)
    {
        double x = model.getNode(i).getX();
        double expected = p0 + (p1 - p0) * x / columnLength;
        double p = nodePressure(model,i);
        QVERIFY2(std::abs(p - expected) <= 1.0e-6 * pScale + 1.0e-9,
                 qPrintable(QString("Node %1 at x = %2: pressure %3 Pa, expected %4 Pa").arg(i).arg(x).arg(p).arg(expected)));
    }

    const RVolume &volume = model.getVolume(0);
    for (uint j=0;j<volume.size();j++)
    {
        RR3Vector v = elementVelocity(model,volume.get(j));
        QVERIFY2(std::abs(v[0] - vx) <= 1.0e-6 * std::abs(vx) &&
                 std::abs(v[1]) <= 1.0e-6 * std::abs(vx) &&
                 std::abs(v[2]) <= 1.0e-6 * std::abs(vx),
                 qPrintable(QString("Element %1: velocity (%2, %3, %4) m/s, expected (%5, 0, 0) m/s")
                            .arg(volume.get(j)).arg(v[0]).arg(v[1]).arg(v[2]).arg(vx)));
    }
}

} // namespace

class TestSolverPotential : public QObject
{
    Q_OBJECT

    private slots:

        //! Pressure difference between the ends drives a uniform Darcy flow.
        void pressureDrivenFlow();

        //! Prescribed inflow velocity gives the matching pressure drop.
        void inflowVelocity();

        //! Prescribed volumetric flow rate is spread over the inlet area.
        void volumetricFlowRate();

        //! Gravity alone drives the flow when both ends are at the same pressure.
        void gravityDrivenFlow();

        //! Closed column under gravity stays at rest with hydrostatic pressure.
        void hydrostaticColumn();

        //! Without a fixed pressure the solution is not unique.
        void missingPressureFails();
};

void TestSolverPotential::pressureDrivenFlow()
{
    const double dp = 1.0e5; // [Pa]

    RModel model = buildColumn();
    setGravity(model,0.0,0.0,0.0);
    setPressure(model,0,dp);
    setPressure(model,1,0.0);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    verifyUniformFlow(model,dp,0.0,conductivity*dp/columnLength);
}

void TestSolverPotential::inflowVelocity()
{
    const double u = 1.0e-4; // [m/s]

    RModel model = buildColumn();
    setGravity(model,0.0,0.0,0.0);

    RBoundaryCondition bc(R_BOUNDARY_CONDITION_INFLOW_VELOCITY);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_VELOCITY)).add(0.0,u);
    model.getSurface(0).addBoundaryCondition(bc);
    setPressure(model,1,0.0);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    verifyUniformFlow(model,u*columnLength/conductivity,0.0,u);
}

void TestSolverPotential::volumetricFlowRate()
{
    const double u = 1.0e-4; // [m/s]
    const double outletPressure = 5.0e3; // [Pa]

    RModel model = buildColumn();
    setGravity(model,0.0,0.0,0.0);

    RBoundaryCondition bc(R_BOUNDARY_CONDITION_INFLOW_VOLURATE);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_VOLUME_FLOW_RATE)).add(0.0,u*columnSide*columnSide);
    model.getSurface(0).addBoundaryCondition(bc);
    setPressure(model,1,outletPressure);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    verifyUniformFlow(model,outletPressure+u*columnLength/conductivity,outletPressure,u);
}

void TestSolverPotential::gravityDrivenFlow()
{
    const double g = 9.80665; // [m/s^2]

    RModel model = buildColumn();
    setGravity(model,-g,0.0,0.0);
    setPressure(model,0,0.0);
    setPressure(model,1,0.0);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    verifyUniformFlow(model,0.0,0.0,-conductivity*density*g);
}

void TestSolverPotential::hydrostaticColumn()
{
    // Vertical column (gravity along -X) closed at the bottom - no flow and
    // hydrostatic pressure p = rho*g*(L - x) below the free top at X = L.
    const double g = 9.80665; // [m/s^2]

    RModel model = buildColumn();
    setGravity(model,-g,0.0,0.0);
    setPressure(model,1,0.0);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double pBottom = density * g * columnLength;
    for (uint i=0;i<model.getNNodes();i++)
    {
        double x = model.getNode(i).getX();
        double expected = density * g * (columnLength - x);
        double p = nodePressure(model,i);
        QVERIFY2(std::abs(p - expected) <= 1.0e-6 * pBottom,
                 qPrintable(QString("Node %1 at x = %2: pressure %3 Pa, expected %4 Pa").arg(i).arg(x).arg(p).arg(expected)));
    }

    const RVolume &volume = model.getVolume(0);
    for (uint j=0;j<volume.size();j++)
    {
        RR3Vector v = elementVelocity(model,volume.get(j));
        QVERIFY2(RR3Vector::euclideanNorm(v) <= 1.0e-6 * conductivity * density * g,
                 qPrintable(QString("Element %1: velocity (%2, %3, %4) m/s, expected none")
                            .arg(volume.get(j)).arg(v[0]).arg(v[1]).arg(v[2])));
    }
}

void TestSolverPotential::missingPressureFails()
{
    RModel model = buildColumn();

    RBoundaryCondition bc(R_BOUNDARY_CONDITION_INFLOW_VELOCITY);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_VELOCITY)).add(0.0,1.0e-4);
    model.getSurface(0).addBoundaryCondition(bc);

    QString solverError;
    QVERIFY(!runSolver(model,solverError));
}

QTEST_APPLESS_MAIN(TestSolverPotential)

#include "tst_solver_potential.moc"

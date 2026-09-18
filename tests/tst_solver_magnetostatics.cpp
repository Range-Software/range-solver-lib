#include <QtTest>

#include <cmath>

#include "rml_model.h"
#include "rsolver.h"
#include "rsolvermagnetostatics.h"

namespace
{

typedef std::array<double,3> Vec;

const double u0 = 1.25663706212e-6;         // [H/m]

constexpr double barLength = 2.0;           // [m]
constexpr double barSide = 0.1;             // [m]
constexpr uint   nSide = 4;                 // divisions across the side
constexpr uint   nSlices = 40;              // divisions along the length
constexpr double currentDensity = 1.0e6;    // [A/m^2]

uint barNodeId(uint i, uint j, uint k)
{
    return (nSide+1)*(nSide+1)*i + (nSide+1)*j + k;
}

uint nBarNodes()
{
    return (nSide+1)*(nSide+1)*(nSlices+1);
}

void addTetrahedron(RModel &model, uint n1, uint n2, uint n3, uint n4)
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
    model.addElement(element,true,0);
}

//! Build a square bar of tetrahedra along the x axis, spanning
//! [0,L] x [0,a] x [0,a], followed by the given free nodes that belong to no
//! element.
RModel buildBar(const QList<Vec> &freeNodes)
{
    RModel model;

    model.setNNodes(nBarNodes() + uint(freeNodes.size()));
    for (uint i=0;i<=nSlices;i++)
    {
        for (uint j=0;j<=nSide;j++)
        {
            for (uint k=0;k<=nSide;k++)
            {
                model.getNode(barNodeId(i,j,k)).set(barLength*double(i)/double(nSlices),
                                                    barSide*double(j)/double(nSide),
                                                    barSide*double(k)/double(nSide));
            }
        }
    }
    for (int i=0;i<freeNodes.size();i++)
    {
        model.getNode(nBarNodes()+uint(i)).set(freeNodes[i][0],freeNodes[i][1],freeNodes[i][2]);
    }

    // Six tetrahedra per brick - a Kuhn decomposition.
    for (uint i=0;i<nSlices;i++)
    {
        for (uint j=0;j<nSide;j++)
        {
            for (uint k=0;k<nSide;k++)
            {
                uint v000 = barNodeId(i,j,k);
                uint v001 = barNodeId(i,j,k+1);
                uint v010 = barNodeId(i,j+1,k);
                uint v011 = barNodeId(i,j+1,k+1);
                uint v100 = barNodeId(i+1,j,k);
                uint v101 = barNodeId(i+1,j,k+1);
                uint v110 = barNodeId(i+1,j+1,k);
                uint v111 = barNodeId(i+1,j+1,k+1);

                addTetrahedron(model,v000,v100,v110,v111);
                addTetrahedron(model,v000,v110,v010,v111);
                addTetrahedron(model,v000,v010,v011,v111);
                addTetrahedron(model,v000,v011,v001,v111);
                addTetrahedron(model,v000,v001,v101,v111);
                addTetrahedron(model,v000,v101,v100,v111);
            }
        }
    }

    model.getVolume(0).setName("Bar");
    model.getTimeSolver().setEnabled(false);

    return model;
}

//! Prescribe a uniform element current density, as the electro-static task
//! would have stored it. The run is marked as a restart, otherwise the solver
//! clears all results - the prescribed current density with them.
void setCurrentDensity(RModel &model, const Vec &J)
{
    uint position = model.addVariable(R_VARIABLE_CURRENT_DENSITY);
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    variable.resize(3,model.getNElements());
    for (uint i=0;i<model.getNElements();i++)
    {
        variable.setValue(0,i,J[0]);
        variable.setValue(1,i,J[1]);
        variable.setValue(2,i,J[2]);
    }
    model.getProblemSetup().setRestart(true);
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

Vec nodeField(const RModel &model, uint nodeID)
{
    uint position = model.findVariable(R_VARIABLE_MAGNETIC_FIELD);
    if (position == RConstants::eod)
    {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan};
    }
    const RVariable &variable = model.getVariable(position);
    return {variable.getValue(0,nodeID),variable.getValue(1,nodeID),variable.getValue(2,nodeID)};
}

double norm(const Vec &v)
{
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double relativeDifference(const Vec &v, const Vec &reference)
{
    Vec d = {v[0]-reference[0],v[1]-reference[1],v[2]-reference[2]};
    return norm(d) / norm(reference);
}

//! Exact field of a straight filament from x = x0 to x = x1 on the line
//! y = yc, z = zc carrying current I in +x.
Vec filamentField(double I, double x0, double x1, double yc, double zc, const Vec &p)
{
    double dy = p[1] - yc;
    double dz = p[2] - zc;
    double d = std::sqrt(dy*dy + dz*dz);
    double c0 = (p[0] - x0) / std::sqrt((p[0]-x0)*(p[0]-x0) + d*d);
    double c1 = (x1 - p[0]) / std::sqrt((x1-p[0])*(x1-p[0]) + d*d);
    double B = u0 * I / (4.0 * RConstants::pi * d) * (c0 + c1);
    // Current along +x, field circulates as x_hat x r_hat.
    return {0.0,-B*dz/d,B*dy/d};
}

//! Gauss-Legendre nodes and weights of order n on [0,1].
void gaussLegendre(uint n, std::vector<double> &x, std::vector<double> &w)
{
    x.resize(n);
    w.resize(n);
    for (uint i=0;i<n;i++)
    {
        // Newton iteration for the roots of the Legendre polynomial on [-1,1].
        double t = std::cos(RConstants::pi * (i + 0.75) / (n + 0.5));
        double dP = 0.0;
        for (uint it=0;it<100;it++)
        {
            double P0 = 1.0;
            double P1 = t;
            for (uint k=2;k<=n;k++)
            {
                double Pk = ((2.0*k-1.0)*t*P1 - (k-1.0)*P0) / k;
                P0 = P1;
                P1 = Pk;
            }
            dP = n * (t*P1 - P0) / (t*t - 1.0);
            double dt = P1 / dP;
            t -= dt;
            if (std::fabs(dt) < 1.0e-15)
            {
                break;
            }
        }
        x[i] = 0.5 * (t + 1.0);
        w[i] = 1.0 / ((1.0 - t*t) * dP * dP);
    }
}

//! Brute force Biot-Savart integral of a uniform current density along +x over
//! the box [0,L] x [y0,y0+w] x [z0,z0+h], midpoint rule on n cells per side.
Vec boxField(double J, double y0, double w, double z0, double h, uint nx, uint ny, uint nz, const Vec &p)
{
    double dx = barLength / nx;
    double dy = w / ny;
    double dz = (nz > 0) ? h / nz : 0.0;
    double dV = dx * dy * ((nz > 0) ? dz : 1.0);
    double By = 0.0;
    double Bz = 0.0;
    for (uint i=0;i<nx;i++)
    {
        double rx = p[0] - (i + 0.5) * dx;
        for (uint j=0;j<ny;j++)
        {
            double ry = p[1] - (y0 + (j + 0.5) * dy);
            for (uint k=0;k<std::max(nz,1u);k++)
            {
                double rz = p[2] - (nz > 0 ? z0 + (k + 0.5) * dz : z0);
                double r2 = rx*rx + ry*ry + rz*rz;
                double f = 1.0 / (r2 * std::sqrt(r2));
                // x_hat x r = (0, -rz, ry)
                By -= rz * f;
                Bz += ry * f;
            }
        }
    }
    double k = u0 * J * dV / (4.0 * RConstants::pi);
    return {0.0,By*k,Bz*k};
}

}

class TestSolverMagnetostatics : public QObject
{
    Q_OBJECT

private slots:

    void segmentMatchesClosedForm();
    void triangleMatchesQuadrature();
    void tetrahedronMatchesQuadrature();
    void barFarFieldMatchesFilament();
    void barNearFieldMatchesDirectIntegral();
    void barFieldVanishesOnAxis();
    void quadratureMatchesExactIntegral();
    void trussMatchesFilament();
    void surfaceStripMatchesDirectIntegral();
    void electrostaticsDrivesMagnetostatics();
    void noCurrentGivesZeroField();

};


void TestSolverMagnetostatics::segmentMatchesClosedForm()
{
    RSolverMagnetostatics::Source source = RSolverMagnetostatics::Source::create(
                2,{Vec{-1.0,0.0,0.0},Vec{1.0,0.0,0.0},Vec{0.0,0.0,0.0},Vec{0.0,0.0,0.0}},{2.0,0.0,0.0});

    Vec p = {0.3,0.5,-0.2};
    Vec B = RSolverMagnetostatics::findSourceField(source,p);
    Vec reference = filamentField(2.0,-1.0,1.0,0.0,0.0,p);
    QVERIFY2(relativeDifference(B,reference) < 1.0e-12,
             qPrintable(QString("Segment field differs by %1").arg(relativeDifference(B,reference))));

    // A point on the line of the segment - inside or beyond it - sees no field.
    for (double x : {-2.0, 0.25, 1.0, 3.0})
    {
        Vec Bl = RSolverMagnetostatics::findSourceField(source,{x,0.0,0.0});
        QCOMPARE(norm(Bl),0.0);
    }
}

void TestSolverMagnetostatics::triangleMatchesQuadrature()
{
    // Sheet current on the triangle (0,0,0), (1,0,0), (0,1,0), against a
    // Gauss-Legendre rule on the triangle collapsed onto the unit square,
    // x = u, y = (1-u)*v.
    const Vec K = {0.3,-1.0,0.0};
    RSolverMagnetostatics::Source source = RSolverMagnetostatics::Source::create(
                3,{Vec{0.0,0.0,0.0},Vec{1.0,0.0,0.0},Vec{0.0,1.0,0.0},Vec{0.0,0.0,0.0}},K);

    std::vector<double> gx, gw;
    gaussLegendre(200,gx,gw);

    auto reference = [&](const Vec &p) -> Vec
    {
        Vec B = {0.0,0.0,0.0};
        for (uint i=0;i<gx.size();i++)
        {
            for (uint j=0;j<gx.size();j++)
            {
                double u = gx[i];
                double v = gx[j];
                Vec r = {p[0]-u,p[1]-(1.0-u)*v,p[2]};
                double r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
                double f = gw[i]*gw[j]*(1.0-u) / (r2*std::sqrt(r2));
                B[0] += (K[1]*r[2]-K[2]*r[1]) * f;
                B[1] += (K[2]*r[0]-K[0]*r[2]) * f;
                B[2] += (K[0]*r[1]-K[1]*r[0]) * f;
            }
        }
        double k = u0 / (4.0 * RConstants::pi);
        return {B[0]*k,B[1]*k,B[2]*k};
    };

    // Above, below, beside and in the plane outside the triangle.
    for (const Vec &p : {Vec{0.2,0.3,0.4},Vec{0.9,0.8,-0.3},Vec{-0.4,0.3,0.2},Vec{1.0,1.0,0.0},Vec{-0.5,-0.2,0.0}})
    {
        Vec B = RSolverMagnetostatics::findSourceFieldExact(source,p);
        Vec Bref = reference(p);
        QVERIFY2(relativeDifference(B,Bref) < 1.0e-6,
                 qPrintable(QString("Triangle field at (%1,%2,%3) differs by %4")
                            .arg(p[0]).arg(p[1]).arg(p[2]).arg(relativeDifference(B,Bref))));
    }

    // Crossing the sheet reverses the tangential field - B = u0*K/2 on either
    // side of an infinite sheet.
    Vec Bup = RSolverMagnetostatics::findSourceFieldExact(source,{0.25,0.25,1.0e-9});
    Vec Bdown = RSolverMagnetostatics::findSourceFieldExact(source,{0.25,0.25,-1.0e-9});
    Vec jump = {Bup[0]-Bdown[0],Bup[1]-Bdown[1],Bup[2]-Bdown[2]};
    Vec expected = {u0*K[1],-u0*K[0],0.0};   // u0 * K x n
    QVERIFY2(relativeDifference(jump,expected) < 1.0e-6,
             qPrintable(QString("Field jump across the sheet differs by %1").arg(relativeDifference(jump,expected))));
}

void TestSolverMagnetostatics::tetrahedronMatchesQuadrature()
{
    const Vec J = {0.3,-1.0,0.5};
    RSolverMagnetostatics::Source source = RSolverMagnetostatics::Source::create(
                4,{Vec{0.0,0.0,0.0},Vec{1.0,0.0,0.0},Vec{0.0,1.0,0.0},Vec{0.0,0.0,1.0}},J);

    QCOMPARE(source.measure,1.0/6.0);

    // Gauss-Legendre rule on the tetrahedron collapsed onto the unit cube,
    // x = u, y = (1-u)*v, z = (1-u)*(1-v)*w.
    std::vector<double> gx, gw;
    gaussLegendre(60,gx,gw);

    auto reference = [&](const Vec &p) -> Vec
    {
        Vec B = {0.0,0.0,0.0};
        for (uint i=0;i<gx.size();i++)
        {
            for (uint j=0;j<gx.size();j++)
            {
                for (uint k=0;k<gx.size();k++)
                {
                    double u = gx[i];
                    double v = gx[j];
                    double w = gx[k];
                    Vec r = {p[0]-u,p[1]-(1.0-u)*v,p[2]-(1.0-u)*(1.0-v)*w};
                    double r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
                    double f = gw[i]*gw[j]*gw[k] * (1.0-u)*(1.0-u)*(1.0-v) / (r2*std::sqrt(r2));
                    B[0] += (J[1]*r[2]-J[2]*r[1]) * f;
                    B[1] += (J[2]*r[0]-J[0]*r[2]) * f;
                    B[2] += (J[0]*r[1]-J[1]*r[0]) * f;
                }
            }
        }
        double k = u0 / (4.0 * RConstants::pi);
        return {B[0]*k,B[1]*k,B[2]*k};
    };

    for (const Vec &p : {Vec{1.0,1.0,1.0},Vec{0.6,0.6,0.6},Vec{-0.3,0.2,0.2},Vec{0.5,-0.5,0.0}})
    {
        Vec B = RSolverMagnetostatics::findSourceFieldExact(source,p);
        Vec Bref = reference(p);
        QVERIFY2(relativeDifference(B,Bref) < 1.0e-6,
                 qPrintable(QString("Tetrahedron field at (%1,%2,%3) differs by %4")
                            .arg(p[0]).arg(p[1]).arg(p[2]).arg(relativeDifference(B,Bref))));
    }

    // The field of a volume current is continuous - at a vertex, on an edge
    // and on a face it has to be the limit of the field next to it.
    for (const Vec &p : {Vec{0.0,0.0,0.0},Vec{0.5,0.0,0.0},Vec{0.2,0.3,0.0},Vec{0.1,0.1,0.1}})
    {
        Vec B = RSolverMagnetostatics::findSourceFieldExact(source,p);
        Vec Bnear = RSolverMagnetostatics::findSourceFieldExact(source,{p[0]-1.0e-7,p[1]-2.0e-7,p[2]-1.5e-7});
        QVERIFY2(std::isfinite(norm(B)) && relativeDifference(B,Bnear) < 1.0e-4,
                 qPrintable(QString("Tetrahedron field at (%1,%2,%3) is not continuous, differs by %4")
                            .arg(p[0]).arg(p[1]).arg(p[2]).arg(relativeDifference(B,Bnear))));
    }

    // Far away the midpoint rule takes over and has to agree.
    Vec pFar = {12.0,-7.0,5.0};
    Vec Bfar = RSolverMagnetostatics::findSourceField(source,pFar);
    Vec BfarExact = RSolverMagnetostatics::findSourceFieldExact(source,pFar);
    QVERIFY2(relativeDifference(Bfar,BfarExact) < 1.0e-3,
             qPrintable(QString("Midpoint rule differs by %1").arg(relativeDifference(Bfar,BfarExact))));
}

void TestSolverMagnetostatics::barFarFieldMatchesFilament()
{
    // Far from a square bar its field is that of a filament on its axis - the
    // first correction of a square cross-section is of fourth order in a/r.
    const double c = barSide / 2.0;
    QList<Vec> points = {
        {barLength/2.0, c + 0.5, c},
        {barLength/2.0, c, c - 0.8},
        {0.3, c - 0.6, c + 0.6}
    };

    RModel model = buildBar(points);
    setCurrentDensity(model,{currentDensity,0.0,0.0});
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double I = currentDensity * barSide * barSide;
    for (int i=0;i<points.size();i++)
    {
        Vec B = nodeField(model,nBarNodes()+uint(i));
        Vec reference = filamentField(I,0.0,barLength,c,c,points[i]);
        QVERIFY2(relativeDifference(B,reference) < 1.0e-3,
                 qPrintable(QString("Far field at point %1 differs by %2").arg(i).arg(relativeDifference(B,reference))));
    }
}

void TestSolverMagnetostatics::barNearFieldMatchesDirectIntegral()
{
    // Close to the bar and on its surface the result has to agree with a direct
    // integral over the conductor.
    const double c = barSide / 2.0;
    QList<Vec> points = {
        {barLength/2.0, barSide + 0.025, c},
        {barLength/2.0 + 0.013, barSide + 0.01, c + 0.02},
        {barLength/2.0, -0.005, -0.005},
        {barLength + 0.5, c + 0.4, c}
    };

    RModel model = buildBar(points);
    setCurrentDensity(model,{currentDensity,0.0,0.0});
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    for (int i=0;i<points.size();i++)
    {
        Vec B = nodeField(model,nBarNodes()+uint(i));
        Vec reference = boxField(currentDensity,0.0,barSide,0.0,barSide,1600,80,80,points[i]);
        QVERIFY2(relativeDifference(B,reference) < 2.0e-3,
                 qPrintable(QString("Near field at point %1 differs by %2").arg(i).arg(relativeDifference(B,reference))));
    }

    // Mesh nodes on the surface - middle of a face and of an edge.
    for (uint nodeID : {barNodeId(nSlices/2,nSide,nSide/2),barNodeId(nSlices/2,0,0)})
    {
        Vec p = {model.getNode(nodeID).getX(),model.getNode(nodeID).getY(),model.getNode(nodeID).getZ()};
        Vec B = nodeField(model,nodeID);
        Vec reference = boxField(currentDensity,0.0,barSide,0.0,barSide,1600,160,160,p);
        QVERIFY2(relativeDifference(B,reference) < 1.0e-2,
                 qPrintable(QString("Surface field at node %1 differs by %2").arg(nodeID).arg(relativeDifference(B,reference))));
    }
}

void TestSolverMagnetostatics::barFieldVanishesOnAxis()
{
    RModel model = buildBar(QList<Vec>());
    setCurrentDensity(model,{currentDensity,0.0,0.0});
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    Vec Baxis = nodeField(model,barNodeId(nSlices/2,nSide/2,nSide/2));
    Vec Bsurface = nodeField(model,barNodeId(nSlices/2,nSide,nSide/2));

    QVERIFY2(norm(Baxis) < 1.0e-6 * norm(Bsurface),
             qPrintable(QString("Field on the axis %1 T against %2 T on the surface").arg(norm(Baxis)).arg(norm(Bsurface))));

    // Every node gets a field, including the ones at the conductor ends.
    const RVariable &variable = model.getVariable(model.findVariable(R_VARIABLE_MAGNETIC_FIELD));
    QCOMPARE(variable.getNValues(),model.getNNodes());
}

void TestSolverMagnetostatics::quadratureMatchesExactIntegral()
{
    // Cells away from the field point are integrated with quadrature rules
    // instead of the closed form - the sum over a whole bar has to stay with
    // the exact one at every node, inside the conductor and on its surface.
    RModel model = buildBar(QList<Vec>());

    // A current density that varies from element to element.
    std::vector<RSolverMagnetostatics::Source> sources;
    for (uint i=0;i<model.getNElements();i++)
    {
        const RElement &element = model.getElement(i);
        std::array<Vec,4> x;
        for (uint k=0;k<4;k++)
        {
            const RNode &node = model.getNode(element.getNodeId(k));
            x[k] = {node.getX(),node.getY(),node.getZ()};
        }
        double s = std::sin(double(i));
        sources.push_back(RSolverMagnetostatics::Source::create(4,x,{currentDensity,0.3*currentDensity*s,-0.2*currentDensity*s*s}));
    }

    double maxError = 0.0;
    double maxField = 0.0;
    for (uint i=0;i<model.getNNodes();i++)
    {
        const RNode &node = model.getNode(i);
        Vec p = {node.getX(),node.getY(),node.getZ()};
        Vec B = {0.0,0.0,0.0};
        Vec Bexact = {0.0,0.0,0.0};
        for (const RSolverMagnetostatics::Source &source : sources)
        {
            Vec dB = RSolverMagnetostatics::findSourceField(source,p);
            Vec dBexact = RSolverMagnetostatics::findSourceFieldExact(source,p);
            for (uint k=0;k<3;k++)
            {
                B[k] += dB[k];
                Bexact[k] += dBexact[k];
            }
        }
        maxError = std::max(maxError,norm({B[0]-Bexact[0],B[1]-Bexact[1],B[2]-Bexact[2]}));
        maxField = std::max(maxField,norm(Bexact));
    }

    QVERIFY2(maxError < 1.0e-4 * maxField,
             qPrintable(QString("Quadrature differs from the exact integral by %1 of the peak field").arg(maxError/maxField)));
}

void TestSolverMagnetostatics::trussMatchesFilament()
{
    const uint nElements = 10;
    const double crossArea = 1.0e-4;

    RModel model;
    model.setNNodes(nElements + 2);
    for (uint i=0;i<=nElements;i++)
    {
        model.getNode(i).set(barLength*double(i)/double(nElements),0.0,0.0);
    }
    Vec p = {0.7,0.2,0.1};
    model.getNode(nElements+1).set(p[0],p[1],p[2]);

    for (uint i=0;i<nElements;i++)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,i);
        element.setNodeId(1,i+1);
        model.addElement(element,true,0);
    }
    model.getLine(0).setCrossArea(crossArea);
    model.getTimeSolver().setEnabled(false);

    setCurrentDensity(model,{currentDensity,0.0,0.0});
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    Vec B = nodeField(model,nElements+1);
    Vec reference = filamentField(currentDensity*crossArea,0.0,barLength,0.0,0.0,p);
    QVERIFY2(relativeDifference(B,reference) < 1.0e-10,
             qPrintable(QString("Truss field differs by %1").arg(relativeDifference(B,reference))));

    // The wire's own nodes lie on its line and see no field from it.
    QCOMPARE(norm(nodeField(model,nElements/2)),0.0);
}

void TestSolverMagnetostatics::surfaceStripMatchesDirectIntegral()
{
    // A thin strip in the plane z = 0, [0,L] x [0,w], meshed with quads.
    const uint nx = 40;
    const uint ny = 4;
    const double width = 0.1;
    const double thickness = 1.0e-3;

    RModel model;
    model.setNNodes((nx+1)*(ny+1) + 1);
    for (uint i=0;i<=nx;i++)
    {
        for (uint j=0;j<=ny;j++)
        {
            model.getNode(i*(ny+1)+j).set(barLength*double(i)/double(nx),width*double(j)/double(ny),0.0);
        }
    }
    Vec p = {barLength/2.0 + 0.01,0.03,0.02};
    model.getNode((nx+1)*(ny+1)).set(p[0],p[1],p[2]);

    for (uint i=0;i<nx;i++)
    {
        for (uint j=0;j<ny;j++)
        {
            RElement element(R_ELEMENT_QUAD1);
            element.setNodeId(0,i*(ny+1)+j);
            element.setNodeId(1,(i+1)*(ny+1)+j);
            element.setNodeId(2,(i+1)*(ny+1)+j+1);
            element.setNodeId(3,i*(ny+1)+j+1);
            model.addElement(element,true,0);
        }
    }
    model.getSurface(0).setThickness(thickness);
    model.getTimeSolver().setEnabled(false);

    setCurrentDensity(model,{currentDensity,0.0,0.0});
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    Vec B = nodeField(model,(nx+1)*(ny+1));
    Vec reference = boxField(currentDensity*thickness,0.0,width,0.0,0.0,4000,200,0,p);
    QVERIFY2(relativeDifference(B,reference) < 2.0e-3,
             qPrintable(QString("Strip field differs by %1").arg(relativeDifference(B,reference))));
}

void TestSolverMagnetostatics::electrostaticsDrivesMagnetostatics()
{
    // A copper bar with a potential difference across its ends - the
    // magneto-static task has to pick up the current of the electro-static one.
    const double conductivity = 5.8e7;  // [S/m]
    const double voltage = 1.0e-3;      // [V]
    const double c = barSide / 2.0;
    Vec p = {barLength/2.0, c + 0.5, c};

    RModel model = buildBar({p});

    RMaterial material(RMaterial::Solid);
    material.setName("Copper");

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,8960.0);
    material.add(densityProperty);

    RMaterialProperty conductivityProperty(RMaterialProperty::ElectricalConductivity);
    conductivityProperty.add(293.15,conductivity);
    material.add(conductivityProperty);

    RMaterialProperty permittivityProperty(RMaterialProperty::RelativePermittivity);
    permittivityProperty.add(293.15,1.0);
    material.add(permittivityProperty);

    model.getVolume(0).setMaterial(material);

    for (uint end=0;end<2;end++)
    {
        uint slice = (end == 0) ? 0 : nSlices;
        for (uint j=0;j<=nSide;j++)
        {
            for (uint k=0;k<=nSide;k++)
            {
                RElement element(R_ELEMENT_POINT);
                element.setNodeId(0,barNodeId(slice,j,k));
                model.addElement(element,true,end);
            }
        }
        RPoint &point = model.getPoint(end);
        point.setVolume(0.0);

        RBoundaryCondition bc(R_BOUNDARY_CONDITION_ELECTRIC_POTENTIAL);
        bc.getComponent(bc.findComponentPosition(R_VARIABLE_ELECTRIC_POTENTIAL)).add(0.0,end == 0 ? voltage : 0.0);
        point.addBoundaryCondition(bc);
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    RProblemTaskItem taskTree;
    taskTree.addChild(RProblemTaskItem(R_PROBLEM_ELECTROSTATICS));
    taskTree.addChild(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));
    model.setProblemTaskTree(taskTree);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double I = conductivity * voltage / barLength * barSide * barSide;
    Vec B = nodeField(model,nBarNodes());
    Vec reference = filamentField(I,0.0,barLength,c,c,p);
    QVERIFY2(relativeDifference(B,reference) < 1.0e-3,
             qPrintable(QString("Field of the electro-static current differs by %1 (B = %2 T, expected %3 T)")
                        .arg(relativeDifference(B,reference)).arg(norm(B)).arg(norm(reference))));
}

void TestSolverMagnetostatics::noCurrentGivesZeroField()
{
    RModel model = buildBar(QList<Vec>());
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_MAGNETOSTATICS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    for (uint i=0;i<model.getNNodes();i++)
    {
        QCOMPARE(norm(nodeField(model,i)),0.0);
    }
}

QTEST_APPLESS_MAIN(TestSolverMagnetostatics)

#include "tst_solver_magnetostatics.moc"

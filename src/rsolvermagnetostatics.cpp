#include <cmath>

#include <omp.h>

#include "rsolvermagnetostatics.h"

namespace
{

typedef std::array<double,3> Vec;

//! A cell closer to the field point than this multiple of its size is
//! integrated exactly, a more distant one with a quadrature rule.
constexpr double nearRatio = 2.0;

//! Beyond this multiple of its size a cell is integrated with the midpoint
//! rule, closer to the field point with a degree two rule.
constexpr double farRatio = 6.0;

inline Vec subtract(const Vec &a, const Vec &b)
{
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

inline Vec cross(const Vec &a, const Vec &b)
{
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

inline double dot(const Vec &a, const Vec &b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double length(const Vec &a)
{
    return std::sqrt(dot(a,a));
}

inline Vec scale(const Vec &a, double f)
{
    return {a[0]*f,a[1]*f,a[2]*f};
}

//! Integrals of a flat triangle (a, b, c) seen from point p.
//!
//! With R the distance from p and m_e the outward in-plane normal of edge e:
//!   I_e   = INT_e 1/R dl              - edge integrals
//!   P_e   = m_e . (x_e - p)           - in-plane distance of p from edge e
//!   d     = n . (p - a)               - height of p above the plane
//!   omega = INT |d|/R^3 dA            - solid angle subtended by the triangle
struct TriangleIntegrals
{
    Vec n;
    Vec m[3];
    double P[3];
    double I[3];
    double d;
    double omega;
};

TriangleIntegrals findTriangleIntegrals(const Vec &a, const Vec &b, const Vec &c, const Vec &p)
{
    TriangleIntegrals ti;

    const Vec *v[3] = {&a,&b,&c};

    Vec n = cross(subtract(b,a),subtract(c,a));
    ti.n = scale(n,1.0/length(n));
    ti.d = dot(ti.n,subtract(p,a));

    for (uint i=0;i<3;i++)
    {
        const Vec &x0 = *v[i];
        const Vec &x1 = *v[(i+1)%3];

        Vec t = subtract(x1,x0);
        double L = length(t);
        t = scale(t,1.0/L);

        // The vertices run counter-clockwise about n, so t x n points out of
        // the triangle.
        ti.m[i] = cross(t,ti.n);
        ti.P[i] = dot(ti.m[i],subtract(x0,p));

        Vec r0 = subtract(x0,p);
        Vec r1 = subtract(x1,p);
        double R0 = length(r0);
        double R1 = length(r1);
        double l0 = dot(r0,t);
        double l1 = dot(r1,t);

        // Squared distance of p from the line of the edge.
        double Rn2 = ti.P[i]*ti.P[i] + ti.d*ti.d;

        if (Rn2 <= 1.0e-24*L*L && l0 <= 0.0 && l1 >= 0.0)
        {
            // p lies on the edge itself and the integral diverges. Callers
            // either multiply it by P = 0 or cancel it against the neighbouring
            // cell that shares the edge.
            ti.I[i] = 0.0;
        }
        else if (l0 + l1 >= 0.0)
        {
            ti.I[i] = std::log((R1 + l1) / (R0 + l0));
        }
        else
        {
            // The same value, written so that nothing cancels behind the edge.
            ti.I[i] = std::log((R0 - l0) / (R1 - l1));
        }
    }

    // Solid angle - Van Oosterom and Strackee.
    Vec r[3] = {subtract(a,p),subtract(b,p),subtract(c,p)};
    double R[3] = {length(r[0]),length(r[1]),length(r[2])};
    double num = dot(r[0],cross(r[1],r[2]));
    double den = R[0]*R[1]*R[2] + dot(r[0],r[1])*R[2] + dot(r[0],r[2])*R[1] + dot(r[1],r[2])*R[0];
    ti.omega = std::fabs(2.0 * std::atan2(num,den));

    return ti;
}

//! G = INT (p - y)/|p - y|^3 dV over a tetrahedron. The integrand is the
//! gradient of 1/|p - y| with respect to y, so the volume integral becomes
//! SUM_faces n_f * INT_f 1/R dA, where
//!   INT_f 1/R dA = SUM_e P_e * I_e - |d| * omega.
Vec findTetrahedronKernel(const std::array<Vec,4> &x, const Vec &p)
{
    // Three vertices of a face followed by the opposite vertex.
    static const uint faces[4][4] = {
        { 1, 2, 3, 0 },
        { 0, 3, 2, 1 },
        { 0, 1, 3, 2 },
        { 0, 2, 1, 3 }
    };

    Vec G = {0.0,0.0,0.0};

    for (const auto &face : faces)
    {
        const Vec &a = x[face[0]];
        Vec b = x[face[1]];
        Vec c = x[face[2]];

        // Orient the face outwards.
        if (dot(cross(subtract(b,a),subtract(c,a)),subtract(x[face[3]],a)) > 0.0)
        {
            std::swap(b,c);
        }

        TriangleIntegrals ti = findTriangleIntegrals(a,b,c,p);

        double phi = -std::fabs(ti.d) * ti.omega;
        for (uint i=0;i<3;i++)
        {
            phi += ti.P[i] * ti.I[i];
        }

        G[0] += ti.n[0] * phi;
        G[1] += ti.n[1] * phi;
        G[2] += ti.n[2] * phi;
    }

    return G;
}

//! G = INT (p - y)/|p - y|^3 dA over a flat triangle. The in-plane part is the
//! surface gradient of 1/|p - y| and becomes SUM_e m_e * I_e, the normal part
//! is the signed solid angle.
Vec findTriangleKernel(const std::array<Vec,4> &x, const Vec &p)
{
    TriangleIntegrals ti = findTriangleIntegrals(x[0],x[1],x[2],p);

    double sign = (ti.d > 0.0) ? 1.0 : ((ti.d < 0.0) ? -1.0 : 0.0);

    Vec G = scale(ti.n,sign*ti.omega);
    for (uint i=0;i<3;i++)
    {
        G[0] += ti.m[i][0] * ti.I[i];
        G[1] += ti.m[i][1] * ti.I[i];
        G[2] += ti.m[i][2] * ti.I[i];
    }

    return G;
}

//! Closed form for a straight segment carrying current I from x0 to x1:
//!   B*4*pi/u0 = I * (a x b) * (|a| + |b|) / (|a|*|b| * (|a|*|b| + a.b))
//! with a = x0 - p and b = x1 - p.
Vec findSegmentKernel(const std::array<Vec,4> &x, const Vec &j, const Vec &p)
{
    Vec B = {0.0,0.0,0.0};

    Vec t = subtract(x[1],x[0]);
    double L = length(t);
    if (L <= 0.0)
    {
        return B;
    }
    double I = dot(j,t) / L;

    Vec a = subtract(x[0],p);
    Vec b = subtract(x[1],p);
    double la = length(a);
    double lb = length(b);
    Vec axb = cross(a,b);
    double lab = la * lb;

    // A point on the line of the segment sees no field from it.
    if (dot(axb,axb) <= 1.0e-24 * lab * lab || lab + dot(a,b) <= 0.0)
    {
        return B;
    }

    return scale(axb,I * (la + lb) / (lab * (lab + dot(a,b))));
}

} // namespace

RSolverMagnetostatics::Source RSolverMagnetostatics::Source::create(uint nVertices, const std::array<std::array<double,3>,4> &x, const std::array<double,3> &j)
{
    Source source;
    source.nVertices = nVertices;
    source.x = x;
    source.j = j;
    source.center = {0.0,0.0,0.0};
    source.size2 = 0.0;
    for (uint i=0;i<nVertices;i++)
    {
        for (uint k=0;k<3;k++)
        {
            source.center[k] += x[i][k] / double(nVertices);
        }
        for (uint l=i+1;l<nVertices;l++)
        {
            Vec e = subtract(x[l],x[i]);
            source.size2 = std::max(source.size2,dot(e,e));
        }
    }
    if (nVertices == 2)
    {
        source.measure = length(subtract(x[1],x[0]));
    }
    else if (nVertices == 3)
    {
        source.measure = 0.5 * length(cross(subtract(x[1],x[0]),subtract(x[2],x[0])));
    }
    else
    {
        source.measure = std::fabs(dot(subtract(x[1],x[0]),cross(subtract(x[2],x[0]),subtract(x[3],x[0])))) / 6.0;
    }
    return source;
}

RSolverMagnetostatics::RSolverMagnetostatics(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
{
    this->problemType = R_PROBLEM_MAGNETOSTATICS;
}

RSolverMagnetostatics::~RSolverMagnetostatics()
{

}

bool RSolverMagnetostatics::hasConverged() const
{
    return true;
}

std::array<double,3> RSolverMagnetostatics::findSourceField(const Source &source, const std::array<double,3> &p)
{
    if (source.nVertices == 2)
    {
        return RSolverMagnetostatics::findSourceFieldExact(source,p);
    }

    Vec r = subtract(p,source.center);
    double r2 = dot(r,r);

    if (r2 < nearRatio*nearRatio*source.size2)
    {
        return RSolverMagnetostatics::findSourceFieldExact(source,p);
    }

    double k = RSolverGeneric::u0 / (4.0 * RConstants::pi);

    if (r2 < farRatio*farRatio*source.size2)
    {
        // Degree two rule - three points on a triangle, four on a tetrahedron.
        const double a = (source.nVertices == 3) ? 2.0/3.0 : 0.5854101966249685;
        const double b = (source.nVertices == 3) ? 1.0/6.0 : 0.1381966011250105;
        const double w = source.measure / double(source.nVertices);

        Vec sum = {0.0,0.0,0.0};
        for (uint l=0;l<source.nVertices;l++)
        {
            sum = {sum[0]+source.x[l][0],sum[1]+source.x[l][1],sum[2]+source.x[l][2]};
        }

        Vec G = {0.0,0.0,0.0};
        for (uint i=0;i<source.nVertices;i++)
        {
            Vec q = {(a-b)*source.x[i][0] + b*sum[0],
                     (a-b)*source.x[i][1] + b*sum[1],
                     (a-b)*source.x[i][2] + b*sum[2]};
            Vec rq = subtract(p,q);
            double rq2 = dot(rq,rq);
            double f = w / (rq2 * std::sqrt(rq2));
            G[0] += rq[0] * f;
            G[1] += rq[1] * f;
            G[2] += rq[2] * f;
        }
        return scale(cross(source.j,G),k);
    }

    // Midpoint rule - B = u0/(4*pi) * J x (p - c) / |p - c|^3 * measure
    return scale(cross(source.j,r),k * source.measure / (r2 * std::sqrt(r2)));
}

std::array<double,3> RSolverMagnetostatics::findSourceFieldExact(const Source &source, const std::array<double,3> &p)
{
    Vec B;

    if (source.nVertices == 2)
    {
        B = findSegmentKernel(source.x,source.j,p);
    }
    else if (source.nVertices == 3)
    {
        B = cross(source.j,findTriangleKernel(source.x,p));
    }
    else
    {
        B = cross(source.j,findTetrahedronKernel(source.x,p));
    }

    return scale(B,RSolverGeneric::u0 / (4.0 * RConstants::pi));
}

void RSolverMagnetostatics::initialize()
{
}

void RSolverMagnetostatics::updateScales()
{

}

void RSolverMagnetostatics::recover()
{
    this->recoverVariable(R_VARIABLE_CURRENT_DENSITY,R_VARIABLE_APPLY_ELEMENT,this->pModel->getNElements(),0,this->elementCurrentDensity.x,0.0);
    this->recoverVariable(R_VARIABLE_CURRENT_DENSITY,R_VARIABLE_APPLY_ELEMENT,this->pModel->getNElements(),1,this->elementCurrentDensity.y,0.0);
    this->recoverVariable(R_VARIABLE_CURRENT_DENSITY,R_VARIABLE_APPLY_ELEMENT,this->pModel->getNElements(),2,this->elementCurrentDensity.z,0.0);
}

void RSolverMagnetostatics::prepare()
{
    this->sources.clear();

    const std::vector<RNode> &nodes = this->pModel->getNodes();

    auto nodePosition = [&nodes](uint nodeID) -> Vec
    {
        return {nodes[nodeID].getX(),nodes[nodeID].getY(),nodes[nodeID].getZ()};
    };

    // An element whose current density is negligible against the largest one
    // is left out - a poorly conducting region such as air would otherwise add
    // a source cell per element and nothing to the field.
    double maxJ2 = 0.0;
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        Vec J = {this->elementCurrentDensity.x[i],this->elementCurrentDensity.y[i],this->elementCurrentDensity.z[i]};
        maxJ2 = std::max(maxJ2,dot(J,J));
    }
    const double minJ2 = 1.0e-20 * maxJ2;

    // Current density of an element scaled by a measure of its entity, or
    // false when the element carries no current.
    auto elementCurrent = [this,minJ2](uint elementID, double factor, Vec &J) -> bool
    {
        J = {this->elementCurrentDensity.x[elementID],
             this->elementCurrentDensity.y[elementID],
             this->elementCurrentDensity.z[elementID]};
        if (dot(J,J) <= minJ2)
        {
            return false;
        }
        J = scale(J,factor);
        return true;
    };

    auto addSource = [this](uint nVertices, const std::array<Vec,4> &x, const Vec &J)
    {
        Source source = Source::create(nVertices,x,J);
        if (source.measure > 0.0)
        {
            this->sources.push_back(source);
        }
    };

    const Vec o = {0.0,0.0,0.0};
    Vec J;

    // Line elements - current J * cross area.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        const RLine &line = this->pModel->getLine(i);
        if (line.getCrossArea() <= 0.0)
        {
            continue;
        }
        for (uint j=0;j<line.size();j++)
        {
            uint elementID = line.get(j);
            const RElement &element = this->pModel->getElement(elementID);
            if (element.getType() != R_ELEMENT_TRUSS1 || !elementCurrent(elementID,line.getCrossArea(),J))
            {
                continue;
            }
            addSource(2,{nodePosition(element.getNodeId(0)),nodePosition(element.getNodeId(1)),o,o},J);
        }
    }

    // Surface elements - sheet current J * thickness. A quadrilateral is split
    // into two triangles.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &surface = this->pModel->getSurface(i);
        if (surface.getThickness() <= 0.0)
        {
            continue;
        }
        for (uint j=0;j<surface.size();j++)
        {
            uint elementID = surface.get(j);
            const RElement &element = this->pModel->getElement(elementID);
            if (!elementCurrent(elementID,surface.getThickness(),J))
            {
                continue;
            }
            if (element.getType() == R_ELEMENT_TRI1)
            {
                addSource(3,{nodePosition(element.getNodeId(0)),
                             nodePosition(element.getNodeId(1)),
                             nodePosition(element.getNodeId(2)),o},J);
            }
            else if (element.getType() == R_ELEMENT_QUAD1)
            {
                addSource(3,{nodePosition(element.getNodeId(0)),
                             nodePosition(element.getNodeId(1)),
                             nodePosition(element.getNodeId(2)),o},J);
                addSource(3,{nodePosition(element.getNodeId(0)),
                             nodePosition(element.getNodeId(2)),
                             nodePosition(element.getNodeId(3)),o},J);
            }
        }
    }

    // Volume elements - current density J.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        const RVolume &volume = this->pModel->getVolume(i);
        for (uint j=0;j<volume.size();j++)
        {
            uint elementID = volume.get(j);
            const RElement &element = this->pModel->getElement(elementID);
            if (element.getType() != R_ELEMENT_TETRA1 || !elementCurrent(elementID,1.0,J))
            {
                continue;
            }
            addSource(4,{nodePosition(element.getNodeId(0)),
                         nodePosition(element.getNodeId(1)),
                         nodePosition(element.getNodeId(2)),
                         nodePosition(element.getNodeId(3))},J);
        }
    }

    if (this->sources.empty())
    {
        RLogger::warning("No element carries a current - the magnetic field is zero. "
                         "The current density is computed by the electro-statics task, which has to run first.\n");
    }
}

void RSolverMagnetostatics::solve()
{
    uint nNodes = this->pModel->getNNodes();

    RLogger::info("Evaluating Biot-Savart integral of %u current carrying cells at %u nodes\n",uint(this->sources.size()),nNodes);

    this->nodeMagneticField.x.resize(nNodes,0.0);
    this->nodeMagneticField.y.resize(nNodes,0.0);
    this->nodeMagneticField.z.resize(nNodes,0.0);

    const std::vector<RNode> &nodes = this->pModel->getNodes();

    // Nodes are processed in blocks so that every source cell is fetched from
    // memory once per block rather than once per node.
    const uint blockSize = 64;
    int64_t nBlocks = (int64_t(nNodes) + blockSize - 1) / blockSize;

    #pragma omp parallel for default(shared) schedule(dynamic,1)
    for (int64_t b=0;b<nBlocks;b++)
    {
        uint first = uint(b) * blockSize;
        uint count = std::min(blockSize,nNodes - first);

        std::array<Vec,blockSize> p;
        std::array<Vec,blockSize> B;
        for (uint i=0;i<count;i++)
        {
            const RNode &node = nodes[first+i];
            p[i] = {node.getX(),node.getY(),node.getZ()};
            B[i] = {0.0,0.0,0.0};
        }

        for (const Source &source : this->sources)
        {
            for (uint i=0;i<count;i++)
            {
                Vec dB = RSolverMagnetostatics::findSourceField(source,p[i]);
                B[i][0] += dB[0];
                B[i][1] += dB[1];
                B[i][2] += dB[2];
            }
        }

        for (uint i=0;i<count;i++)
        {
            this->nodeMagneticField.x[first+i] = B[i][0];
            this->nodeMagneticField.y[first+i] = B[i][1];
            this->nodeMagneticField.z[first+i] = B[i][2];
        }
    }
}

void RSolverMagnetostatics::process()
{

}

void RSolverMagnetostatics::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Magnetic field
    uint magneticFieldPos = this->pModel->findVariable(R_VARIABLE_MAGNETIC_FIELD);
    if (magneticFieldPos == RConstants::eod)
    {
        magneticFieldPos = this->pModel->addVariable(R_VARIABLE_MAGNETIC_FIELD);

        double umin = 0.0;
        double umax = 0.0;
        for (uint i=0;i<this->nodeMagneticField.x.size();i++)
        {
            double u = RR3Vector(this->nodeMagneticField.x[i],
                                 this->nodeMagneticField.y[i],
                                 this->nodeMagneticField.z[i]).length();
            if (i == 0)
            {
                umin = umax = u;
            }
            else
            {
                umin = std::min(umin,u);
                umax = std::max(umax,u);
            }
        }

        this->pModel->getVariable(magneticFieldPos).getVariableData().setMinMaxDisplayValue(umin,umax);
    }
    RVariable &magneticField =  this->pModel->getVariable(magneticFieldPos);

    magneticField.setApplyType(R_VARIABLE_APPLY_NODE);
    magneticField.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        magneticField.setValue(0,i,this->nodeMagneticField.x[i]);
        magneticField.setValue(1,i,this->nodeMagneticField.y[i]);
        magneticField.setValue(2,i,this->nodeMagneticField.z[i]);
    }

    RLogger::unindent();
}

void RSolverMagnetostatics::statistics()
{
    this->printStats(R_VARIABLE_MAGNETIC_FIELD);
    this->processMonitoringPoints();
}

#include <algorithm>
#include <atomic>
#include <omp.h>

#include "rsolverfluid.h"
#include "rmatrixsolver.h"

static const double inv6 = 1.0 / 6.0;

class FluidMatrixContainer
{
    public:

        bool initialized;

        // Element level matricies
        RRMatrix me;    // m
        RRMatrix ce;    // c
        RRMatrix ke;    // k
        RRMatrix ge;    // g
        RRMatrix geT;   // gT
        RRMatrix cpe;   // c+
        RRMatrix cte;   // c~
        RRMatrix ctpe;  // c~+
        RRMatrix kte;   // k~
        RRMatrix ktpe;  // k~+
        RRMatrix ktppe; // k~++
        RRMatrix yte;   // y~
        RRMatrix ytpe;  // y~+
        RRMatrix bte;   // B
        RRMatrix ye;    // y
        RRMatrix ype;   // y+
        RRMatrix the;   // 0
        RRMatrix epe;   // e

        // Element level vectors
        RRVector fv;    // f
        RRVector ftv;   // f~
        RRVector etv;   // e~
        RRVector mv;    // m
        RRVector cv;    // c
        RRVector kv;    // k
        RRVector gv;    // g
        RRVector gvT;   // gT
        RRVector ctv;   // c~
        RRVector ktv;   // k~
        RRVector ytv;   // y~
        RRVector btv;   // B
        RRVector yv;    // y
        RRVector thv;   // 0
        RRVector ev;    // e

        // Element level equations
        RRMatrix Ae11;
        RRMatrix Ae12;
        RRMatrix Ae21;
        RRMatrix Ae22;
        RRVector be1;
        RRVector be2;

        // Assembled element system (thread-local, avoids per-element heap allocation)
        RRMatrix Ae;
        RRVector be;

        // Thread-local temporary vectors (avoid per-element allocation)
        RRVector ax;
        RRVector ay;
        RRVector az;
        RRVector vdiv;

    public:

        FluidMatrixContainer() : initialized(false)
        {

        }

        void resize(uint nen)
        {
            this->me.resize(nen*3,nen*3,0.0);
            this->ce.resize(nen*3,nen*3,0.0);
            this->ke.resize(nen*3,nen*3,0.0);
            this->ge.resize(nen*3,nen,0.0);
            this->geT.resize(nen,nen*3,0.0);
            this->cpe.resize(nen*3,nen*3,0.0);
            this->cte.resize(nen*3,nen*3,0.0);
            this->ctpe.resize(nen*3,nen*3,0.0);
            this->kte.resize(nen*3,nen*3,0.0);
            this->ktpe.resize(nen*3,nen*3,0.0);
            this->ktppe.resize(nen*3,nen*3,0.0);
            this->yte.resize(nen*3,nen,0.0);
            this->ytpe.resize(nen*3,nen*3,0.0);
            this->bte.resize(nen,nen*3,0.0);
            this->ye.resize(nen,nen*3,0.0);
            this->ype.resize(nen,nen*3,0.0);
            this->the.resize(nen,nen,0.0);
            this->epe.resize(nen*3,nen*3,0.0);

            this->fv.resize(nen*3,0.0);
            this->ftv.resize(nen*3,0.0);
            this->etv.resize(nen,0.0);
            this->mv.resize(nen*3,0.0);
            this->cv.resize(nen*3,0.0);
            this->kv.resize(nen*3,0.0);
            this->gv.resize(nen*3,0.0);
            this->gvT.resize(nen,0.0);
            this->ctv.resize(nen*3,0.0);
            this->ktv.resize(nen*3,0.0);
            this->ytv.resize(nen*3,0.0);
            this->btv.resize(nen,0.0);
            this->yv.resize(nen,0.0);
            this->thv.resize(nen,0.0);
            this->ev.resize(nen*3,0.0);

            this->Ae11.resize(nen*3,nen*3,0.0);
            this->Ae12.resize(nen*3,nen,0.0);
            this->Ae21.resize(nen,nen*3,0.0);
            this->Ae22.resize(nen,nen,0.0);
            this->be1.resize(nen*3,0.0);
            this->be2.resize(nen,0.0);

            this->Ae.resize(nen*4,nen*4,0.0);
            this->be.resize(nen*4,0.0);

            // Thread-local temporary vectors
            this->ax.resize(nen,0.0);
            this->ay.resize(nen,0.0);
            this->az.resize(nen,0.0);
            this->vdiv.resize(nen,0.0);

            this->initialized = true;
        }

        void clear()
        {
            this->me.fill(0.0);
            this->ce.fill(0.0);
            this->ke.fill(0.0);
            this->ge.fill(0.0);
            this->geT.fill(0.0);
            this->cpe.fill(0.0);
            this->cte.fill(0.0);
            this->ctpe.fill(0.0);
            this->kte.fill(0.0);
            this->ktpe.fill(0.0);
            this->ktppe.fill(0.0);
            this->yte.fill(0.0);
            this->ytpe.fill(0.0);
            this->bte.fill(0.0);
            this->ye.fill(0.0);
            this->ype.fill(0.0);
            this->the.fill(0.0);
            this->epe.fill(0.0);

            this->fv.fill(0.0);
            this->ftv.fill(0.0);
            this->etv.fill(0.0);
            this->mv.fill(0.0);
            this->cv.fill(0.0);
            this->kv.fill(0.0);
            this->gv.fill(0.0);
            this->gvT.fill(0.0);
            this->ctv.fill(0.0);
            this->ktv.fill(0.0);
            this->ytv.fill(0.0);
            this->btv.fill(0.0);
            this->yv.fill(0.0);
            this->thv.fill(0.0);
            this->ev.fill(0.0);

            this->Ae11.fill(0.0);
            this->Ae12.fill(0.0);
            this->Ae21.fill(0.0);
            this->Ae22.fill(0.0);
            this->be1.fill(0.0);
            this->be2.fill(0.0);

            this->Ae.fill(0.0);
            this->be.fill(0.0);

            // Thread-local temporary vectors
            this->ax.fill(0.0);
            this->ay.fill(0.0);
            this->az.fill(0.0);
            this->vdiv.fill(0.0);
        }
};

bool RSolverFluid::verifyJacobianRequested = false;

const double RSolverFluid::residualDropRatio = 0.1;
// Step of the central difference used by verifyJacobian(). Around the cube root
// of the machine epsilon, which balances the truncation of the difference
// against the cancellation in it. At 1.0e-7 the pressure block of a model at
// rest read 2.5 per cent out, which was the difference and not the matrix.
const double RSolverFluid::differenceStep = 1.0e-5;

const double RSolverFluid::minRelaxation = 0.1;
const double RSolverFluid::relaxationCutFactor = 0.5;
const double RSolverFluid::relaxationGrowFactor = 1.25;
const double RSolverFluid::relaxationRiseTolerance = 0.02;

RSolverFluid::RSolverFluid(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , streamVelocity(1.0)
    , invStreamVelocity(1.0)
    , cvgV(0.0)
    , cvgP(0.0)
    , statsCounter(0)
    , statsOldResidual(0.0)
    , residual(0.0)
    , residualFirst(0.0)
    , previousResidual(0.0)
    , relaxation(1.0)
    , jacobianVerified(false)
    , freezeStabilization(false)
    , xInitialized(false)
{
    this->problemType = R_PROBLEM_FLUID;
    this->nodeAcceleration.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.z.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocity.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocity.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocity.z.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.z.resize(this->pModel->getNNodes(),0.0);
    this->nodePressure.resize(this->pModel->getNNodes(),0.0);
}

RSolverFluid::~RSolverFluid()
{
    this->clearShapeDerivatives();
}

bool RSolverFluid::hasConverged() const
{
    // A task group with no convergence value runs all of its iterations.
    if (this->taskCvgValue <= 0.0)
    {
        return false;
    }

    // Always take a second pass, so that a field which has not moved yet - a
    // model with no inflow, or the very first assembly - is not mistaken for a
    // converged one.
    if (this->taskIteration < 1)
    {
        return false;
    }

    if (this->cvgV >= this->taskCvgValue || this->cvgP >= this->taskCvgValue)
    {
        return false;
    }

    // A small increment on its own is not a solution. A nearly singular system -
    // a model with no pressure reference, or a mesh which can not carry the
    // Reynolds number asked of it - takes tiny steps while its residual stays
    // where it was, or climbs. The residual has to have come down as well.
    if (this->residualFirst > RConstants::eps)
    {
        return (this->residual <= this->residualFirst * RSolverFluid::residualDropRatio);
    }

    // The solve started at a residual of zero - it can only stay there.
    return (this->residual <= RConstants::eps);
}

void RSolverFluid::initialize()
{
    if (!this->xInitialized)
    {
        this->nodeVelocity.x.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocity.y.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocity.z.resize(this->pModel->getNNodes(),0.0);
        this->elementVelocity.x.resize(this->pModel->getNElements(),0.0);
        this->elementVelocity.y.resize(this->pModel->getNElements(),0.0);
        this->elementVelocity.z.resize(this->pModel->getNElements(),0.0);
        this->nodeAcceleration.x.resize(this->pModel->getNNodes(),0.0);
        this->nodeAcceleration.y.resize(this->pModel->getNNodes(),0.0);
        this->nodeAcceleration.z.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocityOld.x.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocityOld.y.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocityOld.z.resize(this->pModel->getNNodes(),0.0);
        this->elementPressure.resize(this->pModel->getNElements(),0.0);
        this->xInitialized = true;
    }
}

void RSolverFluid::updateScales()
{
    RRVector eRo;
    RRVector eU;

    this->generateMaterialVecor(RMaterialProperty::Density,eRo);
    this->generateMaterialVecor(RMaterialProperty::DynamicViscosity,eU);

    this->avgRo = 0.0;
    this->avgU = 0.0;
    uint n = 0;
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (eRo[i] != 0.0 && eU[i] != 0.0)
        {
            this->avgRo += eRo[i];
            this->avgU += eU[i];
            n++;
        }
    }
    if (n > 0)
    {
        this->avgRo /= double(n);
        this->avgU /= double(n);
    }

    this->scales.setMetre(this->findMeshScale());
    this->scales.setSecond(this->findTimeScale());
    this->scales.setKilogram(this->findWeightScale());
}

void RSolverFluid::recover()
{
    this->recoveryStopWatch.reset();
    this->recoveryStopWatch.resume();

    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
//    this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeAcceleration.x);
//    this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeAcceleration.y);
//    this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeAcceleration.z);
    this->recoverVariable(R_VARIABLE_PRESSURE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodePressure,0.0);

    this->recoveryStopWatch.pause();
}

void RSolverFluid::prepare()
{
    RLogger::info("Building matrix system\n");
    RLogger::indent();

    this->buildStopWatch.reset();

    if (this->taskIteration == 0 || this->meshChanged)
    {
        this->generateNodeBook();
        this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);
        this->generateMaterialVecor(RMaterialProperty::DynamicViscosity,this->elementViscosity);

        this->findInputVectors();
    }

    {
        const uint ne = this->pModel->getNElements();
        this->elementPressure.resize(ne,0.0);
        this->elementVelocity.x.resize(ne,0.0);
        this->elementVelocity.y.resize(ne,0.0);
        this->elementVelocity.z.resize(ne,0.0);
        #pragma omp parallel for default(shared)
        for (int64_t i=0;i<int64_t(ne);i++)
        {
            const RElement &el = this->pModel->getElement(uint(i));
            uint nen = el.size();
            double p = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
            for (uint j=0;j<nen;j++)
            {
                uint nid = el.getNodeId(j);
                p  += this->nodePressure[nid];
                vx += this->nodeVelocity.x[nid];
                vy += this->nodeVelocity.y[nid];
                vz += this->nodeVelocity.z[nid];
            }
            double inv = 1.0 / double(nen);
            this->elementPressure[uint(i)]    = p  * inv;
            this->elementVelocity.x[uint(i)] = vx * inv;
            this->elementVelocity.y[uint(i)] = vy * inv;
            this->elementVelocity.z[uint(i)] = vz * inv;
        }
    }

    if (this->meshChanged)
    {
        this->computeElementScales();
        // Node positions changed - cached derivatives are stale.
        this->clearShapeDerivatives();
        this->computeShapeDerivatives();
    }
    if (this->taskIteration == 0)
    {
        this->computeFreePressureNodeHeight();
        this->streamVelocity = RSolverFluid::computeStreamVelocity(*this->pModel,this->nodeVelocity,false);
        this->invStreamVelocity = 1.0 / this->streamVelocity;
    }

    RRVector elementFreePressure;
    RBVector elementFreePressureSetValues;
    this->computeElementFreePressure(elementFreePressure,elementFreePressureSetValues);

    if (this->taskIteration == 0 || this->meshChanged)
    {
        this->buildSparseMatrixPattern(elementFreePressureSetValues);
    }
    this->A.fillValues(0.0);
    this->b.resize(this->nodeBook.getNEnabled());
    this->b.fill(0.0);
    this->x.resize(this->nodeBook.getNEnabled());
    this->x.fill(0.0);

    int np = omp_get_max_threads();

    // Per-thread assembly buffers are members - copying the matrix pattern
    // into them is only needed when the pattern itself was rebuilt.
    std::vector<RSparseMatrix> &Ap = this->threadAssemblyMatrices;
    std::vector<RRVector> &bp = this->threadAssemblyVectors;

    if (this->taskIteration == 0 || this->meshChanged || int(Ap.size()) != np)
    {
        Ap.resize(np);
        bp.resize(np);
        for (int i=0;i<np;i++)
        {
            Ap[i] = this->A;
            bp[i].resize(this->nodeBook.getNEnabled());
        }
    }
    for (int i=0;i<np;i++)
    {
        Ap[i].fillValues(0.0);
        bp[i].fill(0.0);
    }

    std::atomic<bool> abort{false};

    RMatrixManager<FluidMatrixContainer> matrixManager;

    if (this->meshChanged)
    {
        this->elementNormals.resize(this->pModel->getNElements(),RR3Vector(0.0,0.0,0.0));
        this->elementGravityMagnitude.resize(this->pModel->getNElements());
    }

    this->buildStopWatch.resume();

    // Compute element matrices
    #pragma omp parallel for default(shared) private(matrixManager)
    for (int64_t i=0;i<int64_t(this->pModel->getNElements());i++)
    {
        uint elementID = uint(i);

        const RElement &element = this->pModel->getElement(elementID);

        if (abort.load(std::memory_order_relaxed))
        {
            continue;
        }
        try
        {
            uint nInp = RElement::getNIntegrationPoints(element.getType());

            FluidMatrixContainer &mc = matrixManager.getMatricies(element.getType());
            mc.Ae.fill(0.0);
            mc.be.fill(0.0);
            RRMatrix &Ae = mc.Ae;
            RRVector &be = mc.be;

            if (R_ELEMENT_TYPE_IS_SURFACE(element.getType()))
            {
                if (!elementFreePressureSetValues[elementID])
                {
                    continue;
                }

                RR3Vector &normal = this->elementNormals[elementID];
                if (this->meshChanged)
                {
                    // no need to recalculate if mesh does not change !!!
                    element.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]);
                    double gx = elementGravity.x[elementID];
                    double gy = elementGravity.y[elementID];
                    double gz = elementGravity.z[elementID];
                    this->elementGravityMagnitude[elementID] = std::sqrt(gx*gx + gy*gy + gz*gz);
                }

                double ro = this->elementDensity[elementID];

                double fp = elementFreePressure[elementID];
                double gm = this->elementGravityMagnitude[elementID];

                for (uint intPoint=0;intPoint<nInp;intPoint++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),intPoint);
                    const RRVector &N = shapeFunc.getN();
                    double detJ = this->shapeDerivations[elementID]->getJacobian(intPoint);
                    double integValue = detJ * shapeFunc.getW();

                    for (uint m=0;m<element.size();m++)
                    {
                        uint m4 = m*4;

                        double nh = this->freePressureNodeHeight[element.getNodeId(m)];

                        // Pressure vector
                        double value = N[m] * (fp + ro * gm * nh) * integValue;
                        be[m4+0] -= value * normal[0];
                        be[m4+1] -= value * normal[1];
                        be[m4+2] -= value * normal[2];
                    }
                }
                if (this->pModel->getTimeSolver().getEnabled())
                {
                    be *= this->pModel->getTimeSolver().getCurrentTimeStepSize();
                }
            }

            if (R_ELEMENT_TYPE_IS_VOLUME(element.getType()))
            {
                if (!this->computableElements[elementID])
                {
                    continue;
                }
                this->computeElement(elementID,Ae,be,matrixManager);
            }
            this->applyLocalRotations(elementID,Ae);
            this->assemblyMatrix(elementID,Ae,be,Ap[omp_get_thread_num()],bp[omp_get_thread_num()]);
        }
        catch (const RError &rError)
        {
            #pragma omp critical
            {
                RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
            }
            abort.store(true,std::memory_order_relaxed);
        }
    }

#pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->A.getNRows());i++)
    {
        for (int j=0;j<np;j++)
        {
            A.getVector(uint(i)).addVectorValues(Ap[j].getVector(uint(i)));
            this->b[uint(i)] += bp[j][uint(i)];
        }
    }

    this->buildStopWatch.pause();

    if (abort)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
    }

    RLogger::unindent();
}

void RSolverFluid::solve()
{
    RLogger::info("Solving matrix system\n");
    RLogger::indent();

    if (RSolverFluid::verifyJacobianRequested && !this->jacobianVerified)
    {
        this->jacobianVerified = true;
        this->verifyJacobian();
    }

    // The residual assembled by prepare() belongs to the field as it stands
    // now, before this pass moves it, and decides how much of this pass's step
    // is worth taking.
    this->updateResidualAndRelaxation();

    this->solverStopWatch.reset();
    this->solverStopWatch.resume();

    try
    {
        RLogger::indent();
        RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::GMRES));
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,1);
        RLogger::unindent();
    }
    catch (const RError &)
    {
        RLogger::unindent();
        throw;
    }

    this->solverStopWatch.pause();

    this->nodeVelocity.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocity.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocity.z.resize(this->pModel->getNNodes(),0.0);
    this->nodePressure.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.z.resize(this->pModel->getNNodes(),0.0);

    this->updateStopWatch.reset();
    this->updateStopWatch.resume();

    if (!this->pModel->getTimeSolver().getEnabled())
    {
        this->nodeVelocityOld.x = this->nodeVelocity.x;
        this->nodeVelocityOld.y = this->nodeVelocity.y;
        this->nodeVelocityOld.z = this->nodeVelocity.z;
    }

    // Update node velocities and pressures (parallelized).
    // The norm of the increment is accumulated on the way - it is what tells
    // whether the non-linear iteration has settled.
    double dvNorm2 = 0.0;
    double dpNorm2 = 0.0;
    const double omega = this->relaxation;

    #pragma omp parallel for default(shared) reduction(+:dvNorm2,dpNorm2)
    for (int64_t i=0;i<int64_t(this->pModel->getNNodes());i++)
    {
        uint nodeIdx = uint(i);
        uint position = 0;
        double dvx = 0.0;
        double dvy = 0.0;
        double dvz = 0.0;
        double dp = 0.0;

        if (this->nodeBook.getValue(4*nodeIdx+0,position))
        {
            dvx = this->x[position];
        }
        if (this->nodeBook.getValue(4*nodeIdx+1,position))
        {
            dvy = this->x[position];
        }
        if (this->nodeBook.getValue(4*nodeIdx+2,position))
        {
            dvz = this->x[position];
        }
        if (this->nodeBook.getValue(4*nodeIdx+3,position))
        {
            dp = this->x[position];
        }
        if (this->localRotations[nodeIdx].isActive())
        {
            RR3Vector v(dvx,dvy,dvz);
            this->localRotations[nodeIdx].rotateResultsVector(v);
            dvx = v[0];
            dvy = v[1];
            dvz = v[2];
        }

        // Only the relaxed share of the computed step is taken, and it is that
        // share - what the field actually moved by - which the convergence
        // measure below is built from.
        dvx *= omega;
        dvy *= omega;
        dvz *= omega;
        dp  *= omega;

        this->nodeVelocity.x[nodeIdx] += dvx;
        this->nodeVelocity.y[nodeIdx] += dvy;
        this->nodeVelocity.z[nodeIdx] += dvz;
        this->nodePressure[nodeIdx] += dp;

        dvNorm2 += dvx*dvx + dvy*dvy + dvz*dvz;
        dpNorm2 += dp*dp;
    }
    if (this->pModel->getTimeSolver().getEnabled())
    {
        double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
        double invDt = 1.0 / dt;
        #pragma omp parallel for default(shared)
        for (int64_t i=0;i<int64_t(this->pModel->getNNodes());i++)
        {
            uint nodeIdx = uint(i);
            this->nodeAcceleration.x[nodeIdx] = (this->nodeVelocity.x[nodeIdx] - this->nodeVelocityOld.x[nodeIdx]) * invDt;
            this->nodeAcceleration.y[nodeIdx] = (this->nodeVelocity.y[nodeIdx] - this->nodeVelocityOld.y[nodeIdx]) * invDt;
            this->nodeAcceleration.z[nodeIdx] = (this->nodeVelocity.z[nodeIdx] - this->nodeVelocityOld.z[nodeIdx]) * invDt;
        }
    }

    // Compute new velocity norm (parallelized reduction)
    double u = 0.0;
    #pragma omp parallel for default(shared) reduction(+:u)
    for (int64_t i=0;i<int64_t(this->nodeVelocity.x.size());i++)
    {
        double vx = this->nodeVelocity.x[uint(i)];
        double vy = this->nodeVelocity.y[uint(i)];
        double vz = this->nodeVelocity.z[uint(i)];
        u += vx*vx + vy*vy + vz*vz;
    }
    u = std::sqrt(u);
    double p = RRVector::euclideanNorm(this->nodePressure);

    // Relative size of the Newton increment. The field has settled once the
    // step it takes is negligible against the field itself. Both norms are in
    // the same (downscaled) units, so the ratio is dimensionless and needs no
    // scale factor of its own.
    this->cvgV = std::sqrt(dvNorm2) / std::max(u,RConstants::eps);
    this->cvgP = std::sqrt(dpNorm2) / std::max(p,RConstants::eps);

    this->updateStopWatch.pause();

    RLogger::unindent();
}

void RSolverFluid::process()
{

}

void RSolverFluid::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Velocity
    uint velocityPos = this->pModel->findVariable(R_VARIABLE_VELOCITY);
    if (velocityPos == RConstants::eod)
    {
        velocityPos = this->pModel->addVariable(R_VARIABLE_VELOCITY);

        double umin = 0.0;
        double umax = 0.0;
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            double u = RR3Vector(this->nodeVelocity.x[i],
                                 this->nodeVelocity.y[i],
                                 this->nodeVelocity.z[i]).length();
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

        this->pModel->getVariable(velocityPos).getVariableData().setMinMaxDisplayValue(umin,umax);
    }
    RVariable &velocity =  this->pModel->getVariable(velocityPos);

    velocity.setApplyType(R_VARIABLE_APPLY_NODE);
    velocity.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        velocity.setValue(0,i,this->nodeVelocity.x[i]);
        velocity.setValue(1,i,this->nodeVelocity.y[i]);
        velocity.setValue(2,i,this->nodeVelocity.z[i]);
    }

//    // Acceleration
//    uint accelerationPos = this->pModel->findVariable(R_VARIABLE_ACCELERATION);
//    if (accelerationPos == RConstants::eod)
//    {
//        accelerationPos = this->pModel->addVariable(R_VARIABLE_ACCELERATION);

//        double amin = 0.0;
//        double amax = 0.0;
//        for (uint i=0;i<this->nodeVelocity.x.size();i++)
//        {
//            double u = RR3Vector(this->nodeAcceleration.x[i],
//                                 this->nodeAcceleration.y[i],
//                                 this->nodeAcceleration.z[i]).length();
//            if (i == 0)
//            {
//                amin = amax = u;
//            }
//            else
//            {
//                amin = std::min(amin,u);
//                amax = std::max(amax,u);
//            }
//        }

//        this->pModel->getVariable(accelerationPos).getVariableData().setMinMaxDisplayValue(amin,amax);
//    }
//    RVariable &acceleration =  this->pModel->getVariable(accelerationPos);

//    acceleration.setApplyType(R_VARIABLE_APPLY_NODE);
//    acceleration.resize(3,this->pModel->getNNodes());
//    for (uint i=0;i<this->pModel->getNNodes();i++)
//    {
//        acceleration.setValue(0,i,this->nodeAcceleration.x[i]);
//        acceleration.setValue(1,i,this->nodeAcceleration.y[i]);
//        acceleration.setValue(2,i,this->nodeAcceleration.z[i]);
//    }

    // Pressure
    uint pressurePos = this->pModel->findVariable(R_VARIABLE_PRESSURE);
    if (pressurePos == RConstants::eod)
    {
        pressurePos = this->pModel->addVariable(R_VARIABLE_PRESSURE);

        this->pModel->getVariable(pressurePos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodePressure),
                    RStatistics::findMaximumValue(this->nodePressure));
    }
    RVariable &pressure =  this->pModel->getVariable(pressurePos);

    pressure.setApplyType(R_VARIABLE_APPLY_NODE);
    pressure.resize(1,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        pressure.setValue(0,i,this->nodePressure[i]);
    }

    RLogger::unindent();
}

void RSolverFluid::statistics()
{
    // The residual itself was computed at the start of this pass, before the
    // increment was applied - see updateResidualAndRelaxation().
    const double residual = this->residual;
    double convergence = residual - this->statsOldResidual;
    this->statsOldResidual = residual;

    std::vector<RIterationInfoValue> cvgValues;
    cvgValues.push_back(RIterationInfoValue("Solver residual",residual));
    cvgValues.push_back(RIterationInfoValue("Solver convergence",convergence));
    cvgValues.push_back(RIterationInfoValue("Velocity convergence",this->cvgV));
    cvgValues.push_back(RIterationInfoValue("Pressure convergence",this->cvgP));

    RIterationInfo::writeToFile(this->convergenceFileName,this->statsCounter,cvgValues);

    this->printStats(R_VARIABLE_VELOCITY);
    this->printStats(R_VARIABLE_PRESSURE);
    this->processMonitoringPoints();

    RLogger::info("Residual:      % -13g\n",residual);
    RLogger::info("Convergence-R: % -13g\n",convergence);
    RLogger::info("Convergence-V: % -13g\n",this->cvgV);
    RLogger::info("Convergence-P: % -13g\n",this->cvgP);
    if (this->taskCvgValue > 0.0)
    {
        RLogger::info("Relaxation:    % -13g\n",this->relaxation);
        RLogger::info("Residual ratio:% -13g (target %g)\n",
                      (this->residualFirst > RConstants::eps) ? this->residual / this->residualFirst : 0.0,
                      RSolverFluid::residualDropRatio);
        RLogger::info("Convergence target: % -13g%s\n",
                      this->taskCvgValue,
                      this->hasConverged() ? " (reached)" : "");
    }

    RLogger::info("Build time:        %9u [ms]\n",this->buildStopWatch.getMiliSeconds());
    RLogger::info("Solver time:       %9u [ms]\n",this->solverStopWatch.getMiliSeconds());
    RLogger::info("Update time:       %9u [ms]\n",this->updateStopWatch.getMiliSeconds());

    this->statsCounter++;
}

void RSolverFluid::findInputVectors()
{
    RBVector elementVelocitySetValues(this->pModel->getNElements(),false);
    RBVector elementPressureSetValues(this->pModel->getNElements(),false);

    this->elementVelocity.x.resize(this->pModel->getNElements(),0.0);
    this->elementVelocity.y.resize(this->pModel->getNElements(),0.0);
    this->elementVelocity.z.resize(this->pModel->getNElements(),0.0);
    this->nodeAcceleration.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.z.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeVelocityOld.z.resize(this->pModel->getNNodes(),0.0);
    this->elementPressure.resize(this->pModel->getNElements(),0.0);

    // Apply initial conditions
    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        if (this->firstRun)
        {
            // Apply initial conditions
            for (uint j=0;j<pElementGroup->getNInitialConditions();j++)
            {
                RR3Vector velocity(0.0,0.0,0.0);
                double pressure = 0.0;

                bool velocitySet = false;
                bool pressureSet = false;

                const RInitialCondition &ic = pElementGroup->getInitialCondition(j);

                if (ic.getType() == R_INITIAL_CONDITION_VELOCITY)
                {
                    uint icComponentPosition = ic.findComponentPosition(R_VARIABLE_VELOCITY_X);
                    if (icComponentPosition != RConstants::eod)
                    {
                        const RConditionComponent &conditionComponent = ic.getComponent(icComponentPosition);
                        velocity[0] = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                    }

                    icComponentPosition = ic.findComponentPosition(R_VARIABLE_VELOCITY_Y);
                    if (icComponentPosition != RConstants::eod)
                    {
                        const RConditionComponent &conditionComponent = ic.getComponent(icComponentPosition);
                        velocity[1] = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                    }

                    icComponentPosition = ic.findComponentPosition(R_VARIABLE_VELOCITY_Z);
                    if (icComponentPosition != RConstants::eod)
                    {
                        const RConditionComponent &conditionComponent = ic.getComponent(icComponentPosition);
                        velocity[2] = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                    }
                    velocitySet = true;
                }
                else if (ic.getType() == R_INITIAL_CONDITION_PRESSURE)
                {
                    uint icComponentPosition = ic.findComponentPosition(R_VARIABLE_PRESSURE);
                    if (icComponentPosition == RConstants::eod)
                    {
                        continue;
                    }
                    const RConditionComponent &conditionComponent = ic.getComponent(icComponentPosition);
                    pressure = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                    pressureSet = true;
                }

                if (!velocitySet && !pressureSet)
                {
                    continue;
                }

                for (uint k=0;k<pElementGroup->size();k++)
                {
                    if (velocitySet)
                    {
                        this->elementVelocity.x[pElementGroup->get(k)] = velocity[0];
                        this->elementVelocity.y[pElementGroup->get(k)] = velocity[1];
                        this->elementVelocity.z[pElementGroup->get(k)] = velocity[2];
                        elementVelocitySetValues[pElementGroup->get(k)] = true;
                    }
                    if (pressureSet)
                    {
                        this->elementPressure[pElementGroup->get(k)] = pressure;
                        elementPressureSetValues[pElementGroup->get(k)] = true;
                    }
                }
            }
        }
    }

    RBVector elementWall(this->pModel->getNElements(),false);
    RBVector elementFrictionlessWall(this->pModel->getNElements(),false);

    // Apply boundary conditions
    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        // Apply boundary conditions
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            RR3Vector velocity(0.0,0.0,0.0);
            double pressure = 0.0;

            bool wallSet = false;
            bool frictionlessWallSet = false;
            bool velocitySet = false;
            bool pressureSet = false;

            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);

            if (bc.getType() == R_BOUNDARY_CONDITION_WALL)
            {
                wallSet = true;
                velocitySet = true;
            }
            else if (bc.getType() == R_BOUNDARY_CONDITION_WALL_FRICTIONLESS)
            {
                frictionlessWallSet = true;
            }
            else if (bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VELOCITY)
            {
                uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_VELOCITY);
                if (bcComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
                double value = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                RR3Vector normal;

                const RSurface *pSurface = static_cast<const RSurface*>(pElementGroup);
                pSurface->findAverageNormal(this->pModel->getNodes(),this->pModel->getElements(),normal);

                if (pSurface->size() > 0)
                {
                    if (!this->inwardElements[pSurface->get(0)])
                    {
                        normal *= -1.0;
                    }
                }

                velocity[0] = value * normal[0];
                velocity[1] = value * normal[1];
                velocity[2] = value * normal[2];
                velocitySet = true;
            }
            else if (bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VOLURATE)
            {
                uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_VOLUME_FLOW_RATE);
                if (bcComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
                double volumetricFlowRate = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());

                const RSurface *pSurface = static_cast<const RSurface*>(pElementGroup);

                double area = pSurface->findArea(this->pModel->getNodes(),this->pModel->getElements());
                double value = volumetricFlowRate / area;

                RR3Vector normal;
                pSurface->findAverageNormal(this->pModel->getNodes(),this->pModel->getElements(),normal);

                if (pSurface->size() > 0)
                {
                    if (!this->inwardElements[pSurface->get(0)])
                    {
                        normal *= -1.0;
                    }
                }

                velocity[0] = value * normal[0];
                velocity[1] = value * normal[1];
                velocity[2] = value * normal[2];
                velocitySet = true;
            }
            else if (bc.getType() == R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT)
            {
                uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_PRESSURE);
                if (bcComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
                pressure = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
                pressureSet = true;
            }
//            else if (bc.getType() == R_BOUNDARY_CONDITION_PRESSURE_IMPLICIT)
//            {
//                uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_PRESSURE);
//                if (bcComponentPosition == RConstants::eod)
//                {
//                    continue;
//                }
//                const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
//                pressure = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
//                pressureSet = true;
//            }

            if (!wallSet && !frictionlessWallSet && !velocitySet && !pressureSet)
            {
                continue;
            }

            for (uint k=0;k<pElementGroup->size();k++)
            {
                elementWall[pElementGroup->get(k)] = wallSet;
                elementFrictionlessWall[pElementGroup->get(k)] = frictionlessWallSet;

                if (velocitySet)
                {
                    this->elementVelocity.x[pElementGroup->get(k)] = velocity[0];
                    this->elementVelocity.y[pElementGroup->get(k)] = velocity[1];
                    this->elementVelocity.z[pElementGroup->get(k)] = velocity[2];
                    elementVelocitySetValues[pElementGroup->get(k)] = true;
                }
                if (pressureSet)
                {
                    this->elementPressure[pElementGroup->get(k)] = pressure;
                    elementPressureSetValues[pElementGroup->get(k)] = true;
                }
            }
        }
    }

    RBVector elementGravitySetValues;

    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_X,this->elementGravity.x,elementGravitySetValues,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Y,this->elementGravity.y,elementGravitySetValues,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Z,this->elementGravity.z,elementGravitySetValues,true,true,true);

    this->pModel->convertElementToNodeVector(this->elementVelocity.x,elementVelocitySetValues,this->nodeVelocity.x,true);
    this->pModel->convertElementToNodeVector(this->elementVelocity.y,elementVelocitySetValues,this->nodeVelocity.y,true);
    this->pModel->convertElementToNodeVector(this->elementVelocity.z,elementVelocitySetValues,this->nodeVelocity.z,true);
    this->pModel->convertElementToNodeVector(this->elementPressure,elementPressureSetValues,this->nodePressure,true);

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        const RElement &rElement = this->pModel->getElement(i);
        if (elementWall[i])
        {
            for (uint j=0;j<rElement.size();j++)
            {
                this->nodeVelocity.x[rElement.getNodeId(j)] = 0.0;
                this->nodeVelocity.y[rElement.getNodeId(j)] = 0.0;
                this->nodeVelocity.z[rElement.getNodeId(j)] = 0.0;
            }
        }
        if (elementFrictionlessWall[i])
        {
            RR3Vector elementNormal;
            rElement.findNormal(this->pModel->getNodes(),elementNormal[0],elementNormal[1],elementNormal[2]);

            bool hasX = false;
            bool hasY = false;
            bool hasZ = false;

            if (std::fabs(elementNormal[0]) > std::fabs(elementNormal[1]) && std::fabs(elementNormal[0]) > std::fabs(elementNormal[2]))
            {
                hasX = true;
            }
            else
            {
                if (std::fabs(elementNormal[1]) > std::fabs(elementNormal[0]) && std::fabs(elementNormal[1]) > std::fabs(elementNormal[2]))
                {
                    hasY = true;
                }
                else
                {
                    hasZ = true;
                }
            }

            for (uint j=0;j<rElement.size();j++)
            {
                if (hasX)
                {
                    this->nodeVelocity.x[rElement.getNodeId(j)] = 0.0;
                }
                if (hasY)
                {
                    this->nodeVelocity.y[rElement.getNodeId(j)] = 0.0;
                }
                if (hasZ)
                {
                    this->nodeVelocity.z[rElement.getNodeId(j)] = 0.0;
                }
            }
        }
    }

    this->nodeVelocityOld.x = this->nodeVelocity.x;
    this->nodeVelocityOld.y = this->nodeVelocity.y;
    this->nodeVelocityOld.z = this->nodeVelocity.z;

    this->nodeAcceleration.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeAcceleration.z.resize(this->pModel->getNNodes(),0.0);
}

void RSolverFluid::generateNodeBook()
{
    this->nodeBook.resize(this->pModel->getNNodes()*4);
    this->nodeBook.initialize();
    RBVector disabledPositions(this->nodeBook.size(),false);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        bool hasVelocityX = false;
        bool hasVelocityY = false;
        bool hasVelocityZ = false;
        bool hasFriction = false;
        bool hasPressure = false;
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (RBoundaryCondition::getProblemTypeMask(bc.getType()) & R_PROBLEM_FLUID)
            {
                if (bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VELOCITY ||
                    bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VOLURATE ||
                    bc.getType() == R_BOUNDARY_CONDITION_WALL)
                {
                    hasVelocityX = true;
                    hasVelocityY = true;
                    hasVelocityZ = true;
                }
                if (bc.getType() == R_BOUNDARY_CONDITION_WALL_FRICTIONLESS)
                {
                    hasFriction = true;
                }
                if (bc.getType() == R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT)
                {
                    hasPressure = true;
                }
            }
        }
        if (!hasVelocityX && !hasVelocityY && !hasVelocityZ && !hasFriction && !hasPressure)
        {
            continue;
        }
        for (int64_t j=0;j<pElementGroup->size();j++)
        {
            uint elementID = pElementGroup->get(uint(j));
            const RElement &rElement = this->pModel->getElement(elementID);

            bool elementHasVelocityX = hasVelocityX;
            bool elementHasVelocityY = hasVelocityY;
            bool elementHasVelocityZ = hasVelocityZ;

            if (hasFriction)
            {
                RR3Vector elementNormal;
                rElement.findNormal(this->pModel->getNodes(),elementNormal[0],elementNormal[1],elementNormal[2]);

                if (std::fabs(elementNormal[0]) > std::fabs(elementNormal[1]) && std::fabs(elementNormal[0]) > std::fabs(elementNormal[2]))
                {
                    elementHasVelocityX = true;
                }
                else
                {
                    if (std::fabs(elementNormal[1]) > std::fabs(elementNormal[0]) && std::fabs(elementNormal[1]) > std::fabs(elementNormal[2]))
                    {
                        elementHasVelocityY = true;
                    }
                    else
                    {
                        elementHasVelocityZ = true;
                    }
                }
            }
            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                if (elementHasVelocityX)
                {
                    disabledPositions[4*nodeId+0] = true;
                }
                if (elementHasVelocityY)
                {
                    disabledPositions[4*nodeId+1] = true;
                }
                if (elementHasVelocityZ)
                {
                    disabledPositions[4*nodeId+2] = true;
                }
                if (hasPressure)
                {
                    disabledPositions[4*nodeId+3] = true;
                }
            }
        }
    }
    RBVector computableNodes(this->pModel->getNNodes(),false);
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (this->computableElements[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            for (uint j=0;j<rElement.size();j++)
            {
                computableNodes[rElement.getNodeId(j)] = true;
            }
        }
    }
    for (uint i=0;i<computableNodes.size();i++)
    {
        if (!computableNodes[i])
        {
            disabledPositions[4*i+0] = true;
            disabledPositions[4*i+1] = true;
            disabledPositions[4*i+2] = true;
            disabledPositions[4*i+3] = true;
        }
    }

    RSolverGeneric::rebuildNodeBook(this->nodeBook,disabledPositions);
}

void RSolverFluid::computeFreePressureNodeHeight()
{
    this->freePressureNodeHeight.resize(this->pModel->getNNodes());
    this->freePressureNodeHeight.fill(0.0);

    // Based on average gravity vector on surface elements.

    RR3Vector g(0.0,0.0,0.0);

    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);

        if (!rSurface.hasBoundaryCondition(R_BOUNDARY_CONDITION_PRESSURE_IMPLICIT))
        {
            continue;
        }

        if (rSurface.hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION))
        {
            const REnvironmentCondition &ec = rSurface.getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
            uint componentPosition = 0;
            componentPosition = ec.findComponentPosition(R_VARIABLE_G_ACCELERATION_X);
            if (componentPosition != RConstants::eod)
            {
                const RConditionComponent &conditionComponent = ec.getComponent(componentPosition);
                g[0] += conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
            }
            componentPosition = ec.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y);
            if (componentPosition != RConstants::eod)
            {
                const RConditionComponent &conditionComponent = ec.getComponent(componentPosition);
                g[1] += conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
            }
            componentPosition = ec.findComponentPosition(R_VARIABLE_G_ACCELERATION_Z);
            if (componentPosition != RConstants::eod)
            {
                const RConditionComponent &conditionComponent = ec.getComponent(componentPosition);
                g[2] += conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());
            }
        }
    }
    g.normalize();

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        const RNode &rNode = this->pModel->getNode(i);
        this->freePressureNodeHeight[i] = rNode.getX()*g[0]
                                        + rNode.getY()*g[1]
                                        + rNode.getZ()*g[2];
    }

    double md = RStatistics::findMinimumValue(this->freePressureNodeHeight);

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        this->freePressureNodeHeight[i] -= md;
        this->freePressureNodeHeight[i] = std::fabs(this->freePressureNodeHeight[i]);
    }
}

void RSolverFluid::computeShapeDerivatives()
{
    uint ne = this->pModel->getNElements();
    this->shapeDerivations.resize(ne,nullptr);

    #pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(ne);i++)
    {
        uint elementID = uint(i);

        const RElement &rElement = this->pModel->getElement(elementID);
//        if (R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
        {
            if (!this->computableElements[elementID])
            {
                continue;
            }
            if (!this->shapeDerivations[elementID])
            {
                this->shapeDerivations[elementID] = new RElementShapeDerivation(rElement,this->pModel->getNodes(),R_PROBLEM_FLUID);
            }
        }
    }
}

void RSolverFluid::clearShapeDerivatives()
{
    for (uint i=0;i<this->shapeDerivations.size();i++)
    {
        delete this->shapeDerivations[i];
    }
    this->shapeDerivations.clear();
}

void RSolverFluid::computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager)
{
    if (RElement::hasConstantDerivative(this->pModel->getElement(elementID).getType()))
    {
        this->computeElementConstantDerivative(elementID,Ae,be,matrixManager);
    }
    else
    {
        this->computeElementGeneral(elementID,Ae,be,matrixManager);
    }
}

void RSolverFluid::computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();
    uint nInp = RElement::getNIntegrationPoints(element.getType());

    double ro = this->elementDensity[elementID];
    double invro = 1.0 / ro;
    double u = this->elementViscosity[elementID];

    Ae.fill(0.0);
    be.fill(0.0);

    FluidMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());
    matrixCotainer.clear();

    // Element level matricies
    RRMatrix &me = matrixCotainer.me;
    RRMatrix &ce = matrixCotainer.ce;
    RRMatrix &ke = matrixCotainer.ke;
    RRMatrix &ge = matrixCotainer.ge;
    RRMatrix &geT = matrixCotainer.geT;
    RRMatrix &cpe = matrixCotainer.cpe;
    RRMatrix &cte = matrixCotainer.cte;
    RRMatrix &ctpe = matrixCotainer.ctpe;
    RRMatrix &kte = matrixCotainer.kte;
    RRMatrix &ktpe = matrixCotainer.ktpe;
    RRMatrix &ktppe = matrixCotainer.ktppe;
    RRMatrix &yte = matrixCotainer.yte;
    RRMatrix &ytpe = matrixCotainer.ytpe;
    RRMatrix &bte = matrixCotainer.bte;
    RRMatrix &ye = matrixCotainer.ye;
    RRMatrix &ype = matrixCotainer.ype;
    RRMatrix &the = matrixCotainer.the;
    RRMatrix &epe = matrixCotainer.epe;

    // Element level vectors
    RRVector &fv = matrixCotainer.fv;
    RRVector &ftv = matrixCotainer.ftv;
    RRVector &etv = matrixCotainer.etv;
    RRVector &mv = matrixCotainer.mv;
    RRVector &cv = matrixCotainer.cv;
    RRVector &kv = matrixCotainer.kv;
    RRVector &gv = matrixCotainer.gv;
    RRVector &gvT = matrixCotainer.gvT;
    RRVector &ctv = matrixCotainer.ctv;
    RRVector &ktv = matrixCotainer.ktv;
    RRVector &ytv = matrixCotainer.ytv;
    RRVector &btv = matrixCotainer.btv;
    RRVector &yv = matrixCotainer.yv;
    RRVector &thv = matrixCotainer.thv;
    RRVector &ev = matrixCotainer.ev;

    // Element level equations
    RRMatrix &Ae11 = matrixCotainer.Ae11;
    RRMatrix &Ae12 = matrixCotainer.Ae12;
    RRMatrix &Ae21 = matrixCotainer.Ae21;
    RRMatrix &Ae22 = matrixCotainer.Ae22;
    RRVector &be1 = matrixCotainer.be1;
    RRVector &be2 = matrixCotainer.be2;

    double mVh = this->streamVelocity;
    double invmVh = this->invStreamVelocity;
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
    double alphaDt = alpha * dt;

    // Element level input -------------------------------------------
    RR3Vector ve(this->elementVelocity.x[elementID],
                 this->elementVelocity.y[elementID],
                 this->elementVelocity.z[elementID]);
    // Reuse thread-local vectors from matrixContainer
    RRVector &ax = matrixCotainer.ax;
    RRVector &ay = matrixCotainer.ay;
    RRVector &az = matrixCotainer.az;
    for (uint i=0;i<nen;i++)
    {
        ax[i] = this->nodeAcceleration.x[element.getNodeId(i)];
        ay[i] = this->nodeAcceleration.y[element.getNodeId(i)];
        az[i] = this->nodeAcceleration.z[element.getNodeId(i)];
    }
    RR3Vector g(this->elementGravity.x[elementID],
                this->elementGravity.y[elementID],
                this->elementGravity.z[elementID]);
    double p = this->elementPressure[elementID];
    // Element level velocity magnitude and direction, which set the
    // stabilisation parameters and the element length below. Both follow the
    // field being solved. They used to be taken from the velocity of the
    // previous time step, which is the same thing in a steady-state run but in a
    // transient one is the field the step started from, held fixed across every
    // pass of that step - so the stabilisation lagged behind the flow it was
    // meant to stabilise.
    // While the Jacobian check perturbs the field they are held where they
    // were, so that the difference it measures covers only the terms the
    // assembled matrix actually carries - see verifyJacobian().
    double mvh;
    RR3Vector s;
    if (this->freezeStabilization)
    {
        mvh = this->frozenMvh[elementID];
        s = this->frozenS[elementID];
    }
    else
    {
        mvh = ve.length();
        s = ve;
        s.normalize();
    }
    double invmvh = 1.0 / mvh;

    // Reuse thread-local vdiv vector
    RRVector &vdiv = matrixCotainer.vdiv;

    for (uint intPoint=0;intPoint<nInp;intPoint++)
    {
        const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),intPoint);
        const RRVector &N = shapeFunc.getN();
        const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(intPoint);
        double detJ = this->shapeDerivations[elementID]->getJacobian(intPoint);

        double integValue = detJ * shapeFunc.getW();

        // velocity divergence
        vdiv.fill(0.0);
        // partial derivatives
        RR3Vector vex(0.0,0.0,0.0), vey(0.0,0.0,0.0), vez(0.0,0.0,0.0);
        double px = 0.0, py = 0.0, pz = 0.0;
        // element length scale
        double h = 0.0;
        double hn = this->elementScales[elementID];

        for (uint m=0;m<nen;m++)
        {
            vdiv[m] += ve[0] * B[m][0] + ve[1] * B[m][1] + ve[2] * B[m][2];

            vex[0] += B[m][0] * this->nodeVelocity.x[element.getNodeId(m)];
            vex[1] += B[m][0] * this->nodeVelocity.y[element.getNodeId(m)];
            vex[2] += B[m][0] * this->nodeVelocity.z[element.getNodeId(m)];

            vey[0] += B[m][1] * this->nodeVelocity.x[element.getNodeId(m)];
            vey[1] += B[m][1] * this->nodeVelocity.y[element.getNodeId(m)];
            vey[2] += B[m][1] * this->nodeVelocity.z[element.getNodeId(m)];

            vez[0] += B[m][2] * this->nodeVelocity.x[element.getNodeId(m)];
            vez[1] += B[m][2] * this->nodeVelocity.y[element.getNodeId(m)];
            vez[2] += B[m][2] * this->nodeVelocity.z[element.getNodeId(m)];

            px += B[m][0] * this->nodePressure[element.getNodeId(m)];
            py += B[m][1] * this->nodePressure[element.getNodeId(m)];
            pz += B[m][2] * this->nodePressure[element.getNodeId(m)];

            h += std::fabs(s[0]*B[m][0] + s[1]*B[m][1] + s[2]*B[m][2]);
        }
        if (h != 0.0)
        {
            h = 2.0/h;
        }
        // Reynolds numbers
        double roD2u = 0.5 * ro / u;
        double Re = roD2u * mvh * h;
        double Ren = roD2u * mVh * hn;

        // SUPG stabilization parameter
        double Tsupg = 0.0;
        if (mvh > 0.0)
        {
            if (Re > 0.0 && Re <= 3.0)
            {
                Tsupg = h * Re *invmvh * inv6;
            }
            else
            {
                Tsupg = h * invmvh * 0.5;
            }
        }

        // PSPG stabilization parameter
        double Tpspg = 0.0;
        if (mVh > 0.0)
        {
            if (Ren > 0.0 && Ren <= 3.0)
            {
                Tpspg = hn * Ren * invmVh * inv6;
            }
            else
            {
                Tpspg = hn * invmVh * 0.5;
            }
        }

        // LSIC stabilization parameter
        double Tlsic = mvh*h * 0.5;

        double value = 0.0;
        for (uint m=0;m<nen;m++)
        {
            uint64_t m3 = m*3;
            for (uint n=0;n<nen;n++)
            {
                uint64_t n3 = n*3;

                double B00 = B[m][0]*B[n][0];
                double B01 = B[m][0]*B[n][1];
                double B02 = B[m][0]*B[n][2];
                double B10 = B[m][1]*B[n][0];
                double B11 = B[m][1]*B[n][1];
                double B12 = B[m][1]*B[n][2];
                double B20 = B[m][2]*B[n][0];
                double B21 = B[m][2]*B[n][1];
                double B22 = B[m][2]*B[n][2];

                // m matrix
                if (unsteady)
                {
                    value = ro * N[m]*N[n];
                    me[m3+0][m3+0] = value;
                    me[m3+1][m3+1] = value;
                    me[m3+2][m3+2] = value;
                }
                // c matrix
                value = ro * N[m] * vdiv[n];
                ce[m3+0][n3+0] = value;
                ce[m3+1][n3+1] = value;
                ce[m3+2][n3+2] = value;
                // k matrix
                value = ro * u;
                double tmpValue = B00 + B11 + B22;
                ke[m3+0][n3+0] = value * (tmpValue + B00);
                ke[m3+1][n3+0] = value * (           B01);
                ke[m3+2][n3+0] = value * (           B02);
                ke[m3+0][n3+1] = value * (           B10);
                ke[m3+1][n3+1] = value * (tmpValue + B11);
                ke[m3+2][n3+1] = value * (           B12);
                ke[m3+0][n3+2] = value * (           B20);
                ke[m3+1][n3+2] = value * (           B21);
                ke[m3+2][n3+2] = value * (tmpValue + B22);
                // g matrix
                ge[m3+0][n] = B[m][0] * N[n];
                ge[m3+1][n] = B[m][1] * N[n];
                ge[m3+2][n] = B[m][2] * N[n];
                // gT matrix
                geT[n][m3+0] = ge[m3+0][n];
                geT[n][m3+1] = ge[m3+1][n];
                geT[n][m3+2] = ge[m3+2][n];
                // c+ matrix
                value = ro * N[m] * N[n];
                cpe[m3+0][n3+0] = value * vex[0];
                cpe[m3+0][n3+1] = value * vey[0];
                cpe[m3+0][n3+2] = value * vez[0];
                cpe[m3+1][n3+0] = value * vex[1];
                cpe[m3+1][n3+1] = value * vey[1];
                cpe[m3+1][n3+2] = value * vez[1];
                cpe[m3+2][n3+0] = value * vex[2];
                cpe[m3+2][n3+1] = value * vey[2];
                cpe[m3+2][n3+2] = value * vez[2];
                // c~ matrix
                if (unsteady)
                {
                    value = Tsupg * ro * vdiv[m] * N[n];
                    cte[m3+0][n3+0] = value;
                    cte[m3+1][n3+1] = value;
                    cte[m3+2][n3+2] = value;
                }
                // c~+ matrix
                if (unsteady)
                {
                    value = Tsupg * ro * N[n];
                    ctpe[m3+0][n3+0] = value * B[m][0] * ax[n];
                    ctpe[m3+1][n3+0] = value * B[m][0] * ay[n];
                    ctpe[m3+2][n3+0] = value * B[m][0] * az[n];
                    ctpe[m3+0][n3+1] = value * B[m][1] * ax[n];
                    ctpe[m3+1][n3+1] = value * B[m][1] * ay[n];
                    ctpe[m3+2][n3+1] = value * B[m][1] * az[n];
                    ctpe[m3+0][n3+2] = value * B[m][2] * ax[n];
                    ctpe[m3+1][n3+2] = value * B[m][2] * ay[n];
                    ctpe[m3+2][n3+2] = value * B[m][2] * az[n];
                }
                // k~ matrix
                value = Tsupg * ro * vdiv[m] * vdiv[n];
                kte[m3+0][n3+0] = value;
                kte[m3+1][n3+1] = value;
                kte[m3+2][n3+2] = value;
                // k~+ matrix
                value = Tsupg * ro * vdiv[m] * N[n];
                ktpe[m3+0][n3+0] = value * vex[0];
                ktpe[m3+0][n3+1] = value * vey[0];
                ktpe[m3+0][n3+2] = value * vez[0];
                ktpe[m3+1][n3+0] = value * vex[1];
                ktpe[m3+1][n3+1] = value * vey[1];
                ktpe[m3+1][n3+2] = value * vez[1];
                ktpe[m3+2][n3+0] = value * vex[2];
                ktpe[m3+2][n3+1] = value * vey[2];
                ktpe[m3+2][n3+2] = value * vez[2];
                // k~++ matrix
                value = Tsupg * ro * N[n];
                ktppe[m3+0][n3+0] = value * B[m][0] * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0]);
                ktppe[m3+0][n3+1] = value * B[m][1] * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0]);
                ktppe[m3+0][n3+2] = value * B[m][2] * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0]);

                ktppe[m3+1][n3+0] = value * B[m][0] * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1]);
                ktppe[m3+1][n3+1] = value * B[m][1] * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1]);
                ktppe[m3+1][n3+2] = value * B[m][2] * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1]);

                ktppe[m3+2][n3+0] = value * B[m][0] * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]);
                ktppe[m3+2][n3+1] = value * B[m][1] * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]);
                ktppe[m3+2][n3+2] = value * B[m][2] * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]);
                // y~ matrix
                value = Tsupg * vdiv[m];
                yte[m3+0][n] = value * B[n][0];
                yte[m3+1][n] = value * B[n][1];
                yte[m3+2][n] = value * B[n][2];
                // y~+ matrix
                value = Tsupg * N[n];
                ytpe[m3+0][n3+0] = value * B[m][0] * px;
                ytpe[m3+0][n3+1] = value * B[m][1] * px;
                ytpe[m3+0][n3+2] = value * B[m][2] * px;
                ytpe[m3+1][n3+0] = value * B[m][0] * py;
                ytpe[m3+1][n3+1] = value * B[m][1] * py;
                ytpe[m3+1][n3+2] = value * B[m][2] * py;
                ytpe[m3+2][n3+0] = value * B[m][0] * pz;
                ytpe[m3+2][n3+1] = value * B[m][1] * pz;
                ytpe[m3+2][n3+2] = value * B[m][2] * pz;
                // B matrix
                value = Tpspg * N[n];
                bte[m][n3+0] = value * B[m][0];
                bte[m][n3+1] = value * B[m][1];
                bte[m][n3+2] = value * B[m][2];
                // y matrix
                value = Tpspg * vdiv[n];
                ye[m][n3+0] = value * B[m][0];
                ye[m][n3+1] = value * B[m][1];
                ye[m][n3+2] = value * B[m][2];
                // y+ matrix
                value = Tpspg * N[n];
                ype[m][n3+0] = value * (B[m][0]*vex[0] + B[m][1]*vex[1] + B[m][2]*vex[2]);
                ype[m][n3+1] = value * (B[m][0]*vey[0] + B[m][1]*vey[1] + B[m][2]*vey[2]);
                ype[m][n3+2] = value * (B[m][0]*vez[0] + B[m][1]*vez[1] + B[m][2]*vez[2]);
                // 0 matrix
                the[m][n] = Tpspg * ((B00 + B11 + B22) * invro);
                // ep matrix
                value = Tlsic * ro;
                epe[m3+0][n3+0] = B00 * value;
                epe[m3+0][n3+1] = B01 * value;
                epe[m3+0][n3+2] = B02 * value;
                epe[m3+1][n3+0] = B10 * value;
                epe[m3+1][n3+1] = B11 * value;
                epe[m3+1][n3+2] = B12 * value;
                epe[m3+2][n3+0] = B20 * value;
                epe[m3+2][n3+1] = B21 * value;
                epe[m3+2][n3+2] = B22 * value;
            }
            // m vector
            value = N[m] * ro;
            if (unsteady)
            {
                mv[m3+0] = value * ax[m];
                mv[m3+1] = value * ay[m];
                mv[m3+2] = value * az[m];
            }
            // c vector
            cv[m3+0] = value * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0]);
            cv[m3+1] = value * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1]);
            cv[m3+2] = value * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]);
            // k vector
            kv[m3+0] = u * (  B[m][0]*vex[0] + B[m][1]*vey[0] + B[m][2]*vez[0]
                            + B[m][0]*vex[0] + B[m][1]*vex[1] + B[m][2]*vex[2]);
            kv[m3+1] = u * (  B[m][0]*vex[1] + B[m][1]*vey[1] + B[m][2]*vez[1]
                            + B[m][0]*vey[0] + B[m][1]*vey[1] + B[m][2]*vey[2]);
            kv[m3+2] = u * (  B[m][0]*vex[2] + B[m][1]*vey[2] + B[m][2]*vez[2]
                            + B[m][0]*vez[0] + B[m][1]*vez[1] + B[m][2]*vez[2]);
            // g vector
            gv[m3+0] = B[m][0] * p;
            gv[m3+1] = B[m][1] * p;
            gv[m3+2] = B[m][2] * p;
            // gT vector
            gvT[m] = N[m] * (vex[0] + vey[1] + vez[2]);
            // c~ vector
            value = Tsupg * vdiv[m];
            if (unsteady)
            {
                ctv[m3+0] = value * ro * ax[m];
                ctv[m3+1] = value * ro * ay[m];
                ctv[m3+2] = value * ro * az[m];
            }
            // k~ vector
            ktv[m3+0] = value * ro * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0]);
            ktv[m3+1] = value * ro * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1]);
            ktv[m3+2] = value * ro * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]);
            // y~ vector
            ytv[m3+0] = value * px;
            ytv[m3+1] = value * py;
            ytv[m3+2] = value * pz;
            // B vector
            if (unsteady)
            {
                btv[m] = Tpspg * (B[m][0]*ax[m] + B[m][1]*ay[m] + B[m][2]*az[m]);
            }
            // y vector
            yv[m] = Tpspg * (  B[m][0] * (ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0])
                             + B[m][1] * (ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1])
                             + B[m][2] * (ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2]));
            // 0 vector
            thv[m] = Tpspg * ((B[m][0]*px + B[m][1]*py + B[m][2]*pz) * invro);
            // e vector
            value = Tlsic * ro * (vex[0] + vey[1] + vez[2]);
            ev[m3+0] = value * B[m][0];
            ev[m3+1] = value * B[m][1];
            ev[m3+2] = value * B[m][2];
            // f vector
            value = ro * N[m];
            fv[m3+0] = value * g[0];
            fv[m3+1] = value * g[1];
            fv[m3+2] = value * g[2];
            // f~ vector
            value = Tsupg * ro * vdiv[m];
            ftv[m3+0] = value * g[0];
            ftv[m3+1] = value * g[1];
            ftv[m3+2] = value * g[2];
            // e~ vector
            etv[m] = Tpspg * (B[m][0]*g[0] + B[m][1]*g[1] + B[m][2]*g[2]);
        }

        // Assembly element level matrixes
        for (uint m=0;m<nen*3;m++)
        {
            for (uint n=0;n<nen*3;n++)
            {
                if (unsteady)
                {
                    Ae11[m][n] = me[m][n] + cte[m][n]
                               + alphaDt * (
                                              ce[m][n]    + cpe[m][n]
                                            + ctpe[m][n]  + ke[m][n]
                                            + kte[m][n]   + ktpe[m][n]
                                            + ktppe[m][n] + ytpe[m][n] )
                               + dt * epe[m][n];
                }
                else
                {
                    Ae11[m][n] = ce[m][n]    + cpe[m][n]
                               + ke[m][n]
                               + kte[m][n]   + ktpe[m][n]
                               + ktppe[m][n] + ytpe[m][n]
                               + epe[m][n];
                }
            }
            for (uint n=0;n<nen;n++)
            {
                if (unsteady)
                {
                    Ae12[m][n] = -dt * (ge[m][n] + yte[m][n]);
                    Ae21[n][m] = bte[n][m]
                               + dt * geT[n][m]
                               + alphaDt * (ye[n][m] + ype[n][m]);
                }
                else
                {
                    Ae12[m][n] = -ge[m][n] - yte[m][n];
                    Ae21[n][m] =  geT[n][m] + ye[n][m] + ype[n][m];
                }
            }
            if (unsteady)
            {
                be1[m] = dt * (fv[m] + ftv[m] - (  mv[m]  + ctv[m]
                                                 + cv[m]  + kv[m]
                                                 - gv[m]  + ktv[m]
                                                 - ytv[m] + ev[m]) );
            }
            else
            {
                be1[m] = fv[m] + ftv[m] - (  cv[m]  + kv[m]
                                           - gv[m]  + ktv[m]
                                           - ytv[m] + ev[m]);
            }

            for (uint k=0;k<nen;k++)
            {
                for (uint l=0;l<nen;l++)
                {
                    if (unsteady)
                    {
                        Ae22[k][l] = dt * the[k][l];
                    }
                    else
                    {
                        Ae22[k][l] = the[k][l];
                    }
                }
                if (unsteady)
                {
                    be2[k] = dt * (etv[k] - (  btv[k] + gvT[k]
                                             + yv[k] + thv[k]) );
                }
                else
                {
                    be2[k] = etv[k] - (gvT[k] + yv[k] + thv[k]);
                }
            }
        }

        for (uint m=0;m<nen;m++)
        {
            uint64_t m3 = m*3;
            uint64_t m4 = m*4;
            for (uint n=0;n<nen;n++)
            {
                uint64_t n3 = n*3;
                uint64_t n4 = n*4;

                Ae[m4+0][n4+0] += Ae11[m3+0][n3+0] * integValue;
                Ae[m4+0][n4+1] += Ae11[m3+0][n3+1] * integValue;
                Ae[m4+0][n4+2] += Ae11[m3+0][n3+2] * integValue;
                Ae[m4+1][n4+0] += Ae11[m3+1][n3+0] * integValue;
                Ae[m4+1][n4+1] += Ae11[m3+1][n3+1] * integValue;
                Ae[m4+1][n4+2] += Ae11[m3+1][n3+2] * integValue;
                Ae[m4+2][n4+0] += Ae11[m3+2][n3+0] * integValue;
                Ae[m4+2][n4+1] += Ae11[m3+2][n3+1] * integValue;
                Ae[m4+2][n4+2] += Ae11[m3+2][n3+2] * integValue;

                Ae[m4+0][n4+3] += Ae12[m3+0][n] * integValue;
                Ae[m4+1][n4+3] += Ae12[m3+1][n] * integValue;
                Ae[m4+2][n4+3] += Ae12[m3+2][n] * integValue;

                Ae[m4+3][n4+0] += Ae21[m][n3+0] * integValue;
                Ae[m4+3][n4+1] += Ae21[m][n3+1] * integValue;
                Ae[m4+3][n4+2] += Ae21[m][n3+2] * integValue;

                Ae[m4+3][n4+3] += Ae22[m][n] * integValue;
            }
            be[m4+0] += be1[m3+0] * integValue;
            be[m4+1] += be1[m3+1] * integValue;
            be[m4+2] += be1[m3+2] * integValue;

            be[m4+3] += be2[m] * integValue;
        }
    }
}

void RSolverFluid::computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();

    double ro = this->elementDensity[elementID];
    double invro = 1.0 / ro;
    double u = this->elementViscosity[elementID];

    Ae.fill(0.0);
    be.fill(0.0);

    FluidMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());

    double mVh = this->streamVelocity;
    double invmVh = this->invStreamVelocity;
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
    double alphaDt = alpha * dt;

    // Element level input -------------------------------------------
    RR3Vector ve(this->elementVelocity.x[elementID],
                 this->elementVelocity.y[elementID],
                 this->elementVelocity.z[elementID]);
    // Reuse thread-local vectors from matrixContainer
    RRVector &ax = matrixCotainer.ax;
    RRVector &ay = matrixCotainer.ay;
    RRVector &az = matrixCotainer.az;
    for (uint i=0;i<nen;i++)
    {
        ax[i] = this->nodeAcceleration.x[element.getNodeId(i)];
        ay[i] = this->nodeAcceleration.y[element.getNodeId(i)];
        az[i] = this->nodeAcceleration.z[element.getNodeId(i)];
    }

    RR3Vector g(this->elementGravity.x[elementID],
                this->elementGravity.y[elementID],
                this->elementGravity.z[elementID]);
    double p = this->elementPressure[elementID];
    // Element level velocity magnitude and direction, which set the
    // stabilisation parameters and the element length below. Both follow the
    // field being solved. They used to be taken from the velocity of the
    // previous time step, which is the same thing in a steady-state run but in a
    // transient one is the field the step started from, held fixed across every
    // pass of that step - so the stabilisation lagged behind the flow it was
    // meant to stabilise.
    // While the Jacobian check perturbs the field they are held where they
    // were, so that the difference it measures covers only the terms the
    // assembled matrix actually carries - see verifyJacobian().
    double mvh;
    RR3Vector s;
    if (this->freezeStabilization)
    {
        mvh = this->frozenMvh[elementID];
        s = this->frozenS[elementID];
    }
    else
    {
        mvh = ve.length();
        s = ve;
        s.normalize();
    }
    double invmvh = 1.0 / mvh;

    const RRVector &iN = RElement::getMassVector(element.getType());
    const RRMatrix &iNiN = RElement::getMassMatrix(element.getType());
    double wt = RElement::getTotalWeightFactor(element.getType());
    const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(0);

    // Reuse thread-local vdiv vector
    RRVector &vdiv = matrixCotainer.vdiv;
    vdiv.fill(0.0);
    // partial derivatives
    RR3Vector vex(0.0,0.0,0.0), vey(0.0,0.0,0.0), vez(0.0,0.0,0.0);
    double px = 0.0, py = 0.0, pz = 0.0;
    // element length scale
    double h = 0.0;
    double hn = this->elementScales[elementID];

    for (uint m=0;m<nen;m++)
    {
        vdiv[m] += ve[0] * B[m][0] + ve[1] * B[m][1] + ve[2] * B[m][2];

        vex[0] += B[m][0] * this->nodeVelocity.x[element.getNodeId(m)];
        vex[1] += B[m][0] * this->nodeVelocity.y[element.getNodeId(m)];
        vex[2] += B[m][0] * this->nodeVelocity.z[element.getNodeId(m)];

        vey[0] += B[m][1] * this->nodeVelocity.x[element.getNodeId(m)];
        vey[1] += B[m][1] * this->nodeVelocity.y[element.getNodeId(m)];
        vey[2] += B[m][1] * this->nodeVelocity.z[element.getNodeId(m)];

        vez[0] += B[m][2] * this->nodeVelocity.x[element.getNodeId(m)];
        vez[1] += B[m][2] * this->nodeVelocity.y[element.getNodeId(m)];
        vez[2] += B[m][2] * this->nodeVelocity.z[element.getNodeId(m)];

        px += B[m][0] * this->nodePressure[element.getNodeId(m)];
        py += B[m][1] * this->nodePressure[element.getNodeId(m)];
        pz += B[m][2] * this->nodePressure[element.getNodeId(m)];

        h += std::fabs(s[0]*B[m][0] + s[1]*B[m][1] + s[2]*B[m][2]);
    }
    if (h != 0.0)
    {
        h = 2.0/h;
    }
    // Reynolds numbers
    double roD2u = 0.5 * ro / u;
    double Re = roD2u * mvh * h;
    double Ren = roD2u * mVh * hn;

    // SUPG stabilization parameter
    double Tsupg = 0.0;
    if (mvh > 0.0)
    {
        if (Re > 0.0 && Re <= 3.0)
        {
            Tsupg = h * Re * invmvh * inv6;
        }
        else
        {
            Tsupg = h * invmvh * 0.5;
        }
    }

    // PSPG stabilization parameter
    double Tpspg = 0.0;
    if (mVh > 0.0)
    {
        if (Ren > 0.0 && Ren <= 3.0)
        {
            Tpspg = hn * Ren * invmVh * inv6;
        }
        else
        {
            Tpspg = hn * invmVh * 0.5;
        }
    }

    // LSIC stabilization parameter
    double Tlsic = mvh*h * 0.5;

    const double veVex0 = ve[0]*vex[0] + ve[1]*vey[0] + ve[2]*vez[0];
    const double veVex1 = ve[0]*vex[1] + ve[1]*vey[1] + ve[2]*vez[1];
    const double veVex2 = ve[0]*vex[2] + ve[1]*vey[2] + ve[2]*vez[2];
    const double divV = vex[0] + vey[1] + vez[2];

    for (uint m=0;m<nen;m++)
    {
        uint64_t m4 = m*4;
        const double Bm0 = B[m][0];
        const double Bm1 = B[m][1];
        const double Bm2 = B[m][2];
        const double iNm = iN[m];
        const double vdivm = vdiv[m];

        for (uint n=0;n<nen;n++)
        {
            uint64_t n4 = n*4;
            const double Bn0 = B[n][0];
            const double Bn1 = B[n][1];
            const double Bn2 = B[n][2];
            const double iNn = iN[n];
            const double iNiNmn = iNiN[m][n];

            const double B00 = Bm0*Bn0;
            const double B01 = Bm0*Bn1;
            const double B02 = Bm0*Bn2;
            const double B10 = Bm1*Bn0;
            const double B11 = Bm1*Bn1;
            const double B12 = Bm1*Bn2;
            const double B20 = Bm2*Bn0;
            const double B21 = Bm2*Bn1;
            const double B22 = Bm2*Bn2;
            const double tmpValue = B00 + B11 + B22;

            double a00 = 0.0;
            double a01 = 0.0;
            double a02 = 0.0;
            double a10 = 0.0;
            double a11 = 0.0;
            double a12 = 0.0;
            double a20 = 0.0;
            double a21 = 0.0;
            double a22 = 0.0;

            if (unsteady)
            {
                const double mass = ro * iNiNmn;
                const double ctilde = Tsupg * ro * vdivm * iNn;
                const double ctpeScale = Tsupg * ro * iNn;
                a00 += mass + ctilde + alphaDt * (ctpeScale * Bm0 * ax[n]);
                a10 +=                     alphaDt * (ctpeScale * Bm0 * ay[n]);
                a20 +=                     alphaDt * (ctpeScale * Bm0 * az[n]);
                a01 +=                     alphaDt * (ctpeScale * Bm1 * ax[n]);
                a11 += mass + ctilde + alphaDt * (ctpeScale * Bm1 * ay[n]);
                a21 +=                     alphaDt * (ctpeScale * Bm1 * az[n]);
                a02 +=                     alphaDt * (ctpeScale * Bm2 * ax[n]);
                a12 +=                     alphaDt * (ctpeScale * Bm2 * ay[n]);
                a22 += mass + ctilde + alphaDt * (ctpeScale * Bm2 * az[n]);
            }

            const double ce = ro * iNm * vdiv[n];
            a00 += unsteady ? alphaDt * ce : ce;
            a11 += unsteady ? alphaDt * ce : ce;
            a22 += unsteady ? alphaDt * ce : ce;

            const double kScale = ro * u * wt;
            const double k00 = kScale * (tmpValue + B00);
            const double k10 = kScale * B01;
            const double k20 = kScale * B02;
            const double k01 = kScale * B10;
            const double k11 = kScale * (tmpValue + B11);
            const double k21 = kScale * B12;
            const double k02 = kScale * B20;
            const double k12 = kScale * B21;
            const double k22 = kScale * (tmpValue + B22);

            const double cpeScale = ro * iNiNmn;
            const double cpe00 = cpeScale * vex[0];
            const double cpe01 = cpeScale * vey[0];
            const double cpe02 = cpeScale * vez[0];
            const double cpe10 = cpeScale * vex[1];
            const double cpe11 = cpeScale * vey[1];
            const double cpe12 = cpeScale * vez[1];
            const double cpe20 = cpeScale * vex[2];
            const double cpe21 = cpeScale * vey[2];
            const double cpe22 = cpeScale * vez[2];

            const double kte = Tsupg * ro * vdivm * vdiv[n] * wt;
            const double ktpeScale = Tsupg * ro * vdivm * iNn;
            const double ktppeScale = Tsupg * ro * iNn;
            const double ytpeScale = Tsupg * iNn;
            const double epeScale = Tlsic * ro * wt;

            const double steady00 = cpe00 + k00 + kte + ktpeScale * vex[0] + ktppeScale * Bm0 * veVex0 + ytpeScale * Bm0 * px + B00 * epeScale;
            const double steady01 = cpe01 + k01       + ktpeScale * vey[0] + ktppeScale * Bm1 * veVex0 + ytpeScale * Bm1 * px + B01 * epeScale;
            const double steady02 = cpe02 + k02       + ktpeScale * vez[0] + ktppeScale * Bm2 * veVex0 + ytpeScale * Bm2 * px + B02 * epeScale;
            const double steady10 = cpe10 + k10       + ktpeScale * vex[1] + ktppeScale * Bm0 * veVex1 + ytpeScale * Bm0 * py + B10 * epeScale;
            const double steady11 = cpe11 + k11 + kte + ktpeScale * vey[1] + ktppeScale * Bm1 * veVex1 + ytpeScale * Bm1 * py + B11 * epeScale;
            const double steady12 = cpe12 + k12       + ktpeScale * vez[1] + ktppeScale * Bm2 * veVex1 + ytpeScale * Bm2 * py + B12 * epeScale;
            const double steady20 = cpe20 + k20       + ktpeScale * vex[2] + ktppeScale * Bm0 * veVex2 + ytpeScale * Bm0 * pz + B20 * epeScale;
            const double steady21 = cpe21 + k21       + ktpeScale * vey[2] + ktppeScale * Bm1 * veVex2 + ytpeScale * Bm1 * pz + B21 * epeScale;
            const double steady22 = cpe22 + k22 + kte + ktpeScale * vez[2] + ktppeScale * Bm2 * veVex2 + ytpeScale * Bm2 * pz + B22 * epeScale;

            if (unsteady)
            {
                a00 += alphaDt * (steady00 - B00 * epeScale) + dt * B00 * epeScale;
                a01 += alphaDt * (steady01 - B01 * epeScale) + dt * B01 * epeScale;
                a02 += alphaDt * (steady02 - B02 * epeScale) + dt * B02 * epeScale;
                a10 += alphaDt * (steady10 - B10 * epeScale) + dt * B10 * epeScale;
                a11 += alphaDt * (steady11 - B11 * epeScale) + dt * B11 * epeScale;
                a12 += alphaDt * (steady12 - B12 * epeScale) + dt * B12 * epeScale;
                a20 += alphaDt * (steady20 - B20 * epeScale) + dt * B20 * epeScale;
                a21 += alphaDt * (steady21 - B21 * epeScale) + dt * B21 * epeScale;
                a22 += alphaDt * (steady22 - B22 * epeScale) + dt * B22 * epeScale;
            }
            else
            {
                a00 += steady00;
                a01 += steady01;
                a02 += steady02;
                a10 += steady10;
                a11 += steady11;
                a12 += steady12;
                a20 += steady20;
                a21 += steady21;
                a22 += steady22;
            }

            Ae[m4+0][n4+0] += a00;
            Ae[m4+0][n4+1] += a01;
            Ae[m4+0][n4+2] += a02;
            Ae[m4+1][n4+0] += a10;
            Ae[m4+1][n4+1] += a11;
            Ae[m4+1][n4+2] += a12;
            Ae[m4+2][n4+0] += a20;
            Ae[m4+2][n4+1] += a21;
            Ae[m4+2][n4+2] += a22;

            const double ge0 = Bm0 * iNn;
            const double ge1 = Bm1 * iNn;
            const double ge2 = Bm2 * iNn;
            const double yteScale = Tsupg * vdivm * wt;
            const double ae12Scale = unsteady ? -dt : -1.0;
            Ae[m4+0][n4+3] += ae12Scale * (ge0 + yteScale * Bn0);
            Ae[m4+1][n4+3] += ae12Scale * (ge1 + yteScale * Bn1);
            Ae[m4+2][n4+3] += ae12Scale * (ge2 + yteScale * Bn2);

            const double bteScale = Tpspg * iNn;
            const double geTScale = iNm;
            const double yeScale = Tpspg * vdiv[n] * wt;
            const double ypeScale = Tpspg * iNn;
            const double ype0 = ypeScale * (Bm0*vex[0] + Bm1*vex[1] + Bm2*vex[2]);
            const double ype1 = ypeScale * (Bm0*vey[0] + Bm1*vey[1] + Bm2*vey[2]);
            const double ype2 = ypeScale * (Bm0*vez[0] + Bm1*vez[1] + Bm2*vez[2]);
            if (unsteady)
            {
                Ae[m4+3][n4+0] += bteScale * Bm0 + dt * (Bn0 * geTScale) + alphaDt * (yeScale * Bm0 + ype0);
                Ae[m4+3][n4+1] += bteScale * Bm1 + dt * (Bn1 * geTScale) + alphaDt * (yeScale * Bm1 + ype1);
                Ae[m4+3][n4+2] += bteScale * Bm2 + dt * (Bn2 * geTScale) + alphaDt * (yeScale * Bm2 + ype2);
            }
            else
            {
                Ae[m4+3][n4+0] += Bn0 * geTScale + yeScale * Bm0 + ype0;
                Ae[m4+3][n4+1] += Bn1 * geTScale + yeScale * Bm1 + ype1;
                Ae[m4+3][n4+2] += Bn2 * geTScale + yeScale * Bm2 + ype2;
            }

            const double theta = Tpspg * (tmpValue * invro) * wt;
            Ae[m4+3][n4+3] += unsteady ? dt * theta : theta;
        }

        const double cv0 = iNm * ro * veVex0;
        const double cv1 = iNm * ro * veVex1;
        const double cv2 = iNm * ro * veVex2;
        const double kv0 = u * (Bm0*vex[0] + Bm1*vey[0] + Bm2*vez[0] + Bm0*vex[0] + Bm1*vex[1] + Bm2*vex[2]) * wt;
        const double kv1 = u * (Bm0*vex[1] + Bm1*vey[1] + Bm2*vez[1] + Bm0*vey[0] + Bm1*vey[1] + Bm2*vey[2]) * wt;
        const double kv2 = u * (Bm0*vex[2] + Bm1*vey[2] + Bm2*vez[2] + Bm0*vez[0] + Bm1*vez[1] + Bm2*vez[2]) * wt;
        const double gv0 = Bm0 * p * wt;
        const double gv1 = Bm1 * p * wt;
        const double gv2 = Bm2 * p * wt;
        const double ktvScale = Tsupg * ro * vdivm * wt;
        const double ktv0 = ktvScale * veVex0;
        const double ktv1 = ktvScale * veVex1;
        const double ktv2 = ktvScale * veVex2;
        const double ytvScale = Tsupg * vdivm * wt;
        const double ytv0 = ytvScale * px;
        const double ytv1 = ytvScale * py;
        const double ytv2 = ytvScale * pz;
        const double evScale = Tlsic * ro * divV * wt;
        const double ev0 = evScale * Bm0;
        const double ev1 = evScale * Bm1;
        const double ev2 = evScale * Bm2;
        const double fvScale = ro * iNm;
        const double ftvScale = Tsupg * ro * vdivm * wt;

        if (unsteady)
        {
            const double mvScale = iNm * ro;
            const double ctvScale = Tsupg * ro * vdivm * wt;
            be[m4+0] += dt * (fvScale * g[0] + ftvScale * g[0] - (mvScale * ax[m] + ctvScale * ax[m] + cv0 + kv0 - gv0 + ktv0 - ytv0 + ev0));
            be[m4+1] += dt * (fvScale * g[1] + ftvScale * g[1] - (mvScale * ay[m] + ctvScale * ay[m] + cv1 + kv1 - gv1 + ktv1 - ytv1 + ev1));
            be[m4+2] += dt * (fvScale * g[2] + ftvScale * g[2] - (mvScale * az[m] + ctvScale * az[m] + cv2 + kv2 - gv2 + ktv2 - ytv2 + ev2));
            be[m4+3] += dt * (Tpspg * (Bm0*g[0] + Bm1*g[1] + Bm2*g[2]) * wt
                             - (Tpspg * (Bm0*ax[m] + Bm1*ay[m] + Bm2*az[m]) * wt
                                + iNm * divV
                                + Tpspg * (Bm0*veVex0 + Bm1*veVex1 + Bm2*veVex2) * wt
                                + Tpspg * ((Bm0*px + Bm1*py + Bm2*pz) * invro) * wt));
        }
        else
        {
            be[m4+0] += fvScale * g[0] + ftvScale * g[0] - (cv0 + kv0 - gv0 + ktv0 - ytv0 + ev0);
            be[m4+1] += fvScale * g[1] + ftvScale * g[1] - (cv1 + kv1 - gv1 + ktv1 - ytv1 + ev1);
            be[m4+2] += fvScale * g[2] + ftvScale * g[2] - (cv2 + kv2 - gv2 + ktv2 - ytv2 + ev2);
            be[m4+3] += Tpspg * (Bm0*g[0] + Bm1*g[1] + Bm2*g[2]) * wt
                       - (iNm * divV
                          + Tpspg * (Bm0*veVex0 + Bm1*veVex1 + Bm2*veVex2) * wt
                          + Tpspg * ((Bm0*px + Bm1*py + Bm2*pz) * invro) * wt);
        }
    }

    double detJ = this->shapeDerivations[elementID]->getJacobian(0);
    Ae *= detJ;
    be *= detJ;
}

void RSolverFluid::setVerifyJacobian(bool verifyJacobian)
{
    RSolverFluid::verifyJacobianRequested = verifyJacobian;
}

void RSolverFluid::updateNodeAcceleration()
{
    if (!this->pModel->getTimeSolver().getEnabled())
    {
        return;
    }

    const double invDt = 1.0 / this->pModel->getTimeSolver().getCurrentTimeStepSize();

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        this->nodeAcceleration.x[i] = (this->nodeVelocity.x[i] - this->nodeVelocityOld.x[i]) * invDt;
        this->nodeAcceleration.y[i] = (this->nodeVelocity.y[i] - this->nodeVelocityOld.y[i]) * invDt;
        this->nodeAcceleration.z[i] = (this->nodeVelocity.z[i] - this->nodeVelocityOld.z[i]) * invDt;
    }
}

void RSolverFluid::verifyJacobian()
{
    const uint nEnabled = this->nodeBook.getNEnabled();
    const uint nNodes = this->pModel->getNNodes();
    const uint nElements = this->pModel->getNElements();

    // Every reassembly below has to see the field as this routine leaves it.
    // On the first pass of a solve, and whenever the mesh has just changed,
    // prepare() rebuilds the field from the boundary conditions and the node
    // book along with it, which would wipe out each perturbation and could even
    // change the size of the system in the middle of the sweep. Both are held
    // for the duration - neither the mesh nor the conditions move here.
    const uint taskIterationBackup = this->taskIteration;
    const bool meshChangedBackup = this->meshChanged;
    this->taskIteration = std::max(this->taskIteration,uint(1));
    this->meshChanged = false;

    RLogger::info("Verifying the Jacobian against a finite difference of the residual\n");
    RLogger::indent();
    RLogger::info("Unknowns:      %u\n",nEnabled);
    RLogger::info("Assemblies:    %u\n",2*nEnabled);
    if (nElements > 0)
    {
        RLogger::info("Density:       %g\n",this->elementDensity[0]);
        RLogger::info("Viscosity:     %g\n",this->elementViscosity[0]);
        RLogger::info("Density*visc.: %g\n",this->elementDensity[0]*this->elementViscosity[0]);
    }

    // Map every enabled row of the matrix back to the node and the component
    // it stands for - 0,1,2 are the velocity components, 3 is the pressure.
    std::vector<uint> dofNode(nEnabled,0);
    std::vector<uint> dofComponent(nEnabled,0);
    for (uint i=0;i<nNodes;i++)
    {
        for (uint c=0;c<4;c++)
        {
            uint position = 0;
            if (this->nodeBook.getValue(4*i+c,position))
            {
                dofNode[position] = i;
                dofComponent[position] = c;
            }
        }
    }

    // The system as it was assembled from the unperturbed field.
    const RSparseMatrix Aref(this->A);
    const RRVector bRef(this->b);

    // Hold the stabilisation parameters at the values they have for that same
    // field. They follow the velocity, but the matrix carries no derivative of
    // them, so letting them move would put an expected disagreement into every
    // term they touch and drown the one being looked for.
    this->frozenMvh.resize(nElements,0.0);
    this->frozenS.resize(nElements,RR3Vector(0.0,0.0,0.0));
    for (uint i=0;i<nElements;i++)
    {
        RR3Vector ve(this->elementVelocity.x[i],this->elementVelocity.y[i],this->elementVelocity.z[i]);
        this->frozenMvh[i] = ve.length();
        ve.normalize();
        this->frozenS[i] = ve;
    }
    this->freezeStabilization = true;

    // Ratios of the assembled entry to the measured one, kept per block of the
    // system so that a factor sitting on one block shows up as a factor.
    const uint nBlocks = 4;
    const QString blockName[nBlocks] = { "velocity / velocity",
                                         "velocity / pressure",
                                         "pressure / velocity",
                                         "pressure / pressure" };
    std::vector<std::vector<double>> blockRatios(nBlocks);

    double worstDifference = 0.0;
    uint worstRow = 0;
    uint worstColumn = 0;
    double worstAssembled = 0.0;
    double worstMeasured = 0.0;

    RRVector bPlus;
    RRVector bMinus;

    for (uint j=0;j<nEnabled;j++)
    {
        const uint node = dofNode[j];
        const uint component = dofComponent[j];

        double *field = nullptr;
        switch (component)
        {
            case 0:  field = &this->nodeVelocity.x[node]; break;
            case 1:  field = &this->nodeVelocity.y[node]; break;
            case 2:  field = &this->nodeVelocity.z[node]; break;
            default: field = &this->nodePressure[node];   break;
        }

        const double value = *field;
        const double step = RSolverFluid::differenceStep * std::max(std::fabs(value),1.0);

        *field = value + step;
        this->updateNodeAcceleration();
        this->prepare();
        bPlus = this->b;

        *field = value - step;
        this->updateNodeAcceleration();
        this->prepare();
        bMinus = this->b;

        *field = value;
        this->updateNodeAcceleration();

        for (uint i=0;i<nEnabled;i++)
        {
            // The iteration solves A*dx = b and applies x += dx, so it drives
            // b(x) to zero and the exact Newton matrix is -db/dx.
            const double measured = -(bPlus[i] - bMinus[i]) / (2.0 * step);
            const double assembled = Aref.findValue(i,j);

            const double scale = std::max(std::fabs(measured),std::fabs(assembled));
            if (scale < RConstants::eps)
            {
                continue;
            }

            const double difference = std::fabs(assembled - measured) / scale;
            if (difference > worstDifference)
            {
                worstDifference = difference;
                worstRow = i;
                worstColumn = j;
                worstAssembled = assembled;
                worstMeasured = measured;
            }

            if (std::fabs(measured) < RConstants::eps)
            {
                continue;
            }

            const uint block = (dofComponent[i] < 3 ? 0 : 2) + (component < 3 ? 0 : 1);
            blockRatios[block].push_back(assembled / measured);
        }
    }

    this->freezeStabilization = false;
    this->taskIteration = taskIterationBackup;
    this->meshChanged = meshChangedBackup;
    this->A = Aref;
    this->b = bRef;

    RLogger::info("Ratio of the assembled entry to the measured one, by block:\n");
    RLogger::indent();
    for (uint i=0;i<nBlocks;i++)
    {
        std::vector<double> &ratios = blockRatios[i];
        if (ratios.empty())
        {
            RLogger::info("%-20s no entries\n",blockName[i].toUtf8().constData());
            continue;
        }
        std::sort(ratios.begin(),ratios.end());
        RLogger::info("%-20s entries %6u | min % -12g | median % -12g | max % -12g\n",
                      blockName[i].toUtf8().constData(),
                      uint(ratios.size()),
                      ratios.front(),
                      ratios[ratios.size()/2],
                      ratios.back());

        // Entries which agree sit on one value, so grouping what is within the
        // noise of the difference tells the terms of the block apart: one group
        // at 1 and one at something else names a factor and says how much of
        // the block carries it.
        std::vector<double> groupValue;
        std::vector<uint> groupCount;
        for (uint j=0;j<ratios.size();j++)
        {
            if (!groupValue.empty() &&
                std::fabs(ratios[j] - groupValue.back()) <= 0.02 * std::max(std::fabs(groupValue.back()),1.0e-3))
            {
                groupCount.back()++;
                continue;
            }
            groupValue.push_back(ratios[j]);
            groupCount.push_back(1);
        }

        RLogger::indent();
        for (uint j=0;j<groupValue.size() && j<6;j++)
        {
            RLogger::info("%6u entries near % -12g\n",groupCount[j],groupValue[j]);
        }
        if (groupValue.size() > 6)
        {
            RLogger::info("%6u further groups\n",uint(groupValue.size())-6);
        }
        RLogger::unindent(false);
    }
    RLogger::unindent(false);

    RLogger::info("Largest relative disagreement: %g\n",worstDifference);
    RLogger::indent();
    RLogger::info("row    node %u component %u\n",dofNode[worstRow],dofComponent[worstRow]);
    RLogger::info("column node %u component %u\n",dofNode[worstColumn],dofComponent[worstColumn]);
    RLogger::info("assembled % -13g measured % -13g\n",worstAssembled,worstMeasured);
    RLogger::unindent(false);

    RLogger::info("A block whose ratio is 1 throughout is assembled correctly. A block\n");
    RLogger::info("whose ratio is the same number other than 1 throughout carries that\n");
    RLogger::info("factor too many. A scattered ratio means the term is wrong in form.\n");

    RLogger::unindent();
}

void RSolverFluid::updateResidualAndRelaxation()
{
    const double secondScale = this->scales.getSecond();
    const double scale = (secondScale * secondScale) / this->scales.getKilogram();

    this->previousResidual = this->residual;
    this->residual = RRVector::euclideanNorm(this->b) * scale;

    if (this->taskIteration == 0)
    {
        // Each solve - each time step of a transient run - starts its own
        // history and takes the full step until told otherwise.
        this->residualFirst = this->residual;
        this->previousResidual = this->residual;
        this->relaxation = 1.0;
        return;
    }

    // The residual of this pass judges the step the previous pass took: the
    // matrix is not the exact derivative of the residual, so a full step can
    // overshoot and leave the field further from a solution than it started.
    // Retreat quickly when that happens and return to the full step slowly.
    // Only a rise worth reacting to counts. The residual of a stabilised flow
    // model wanders a little from pass to pass even while it is converging, and
    // retreating on every one of those stalls the iteration at a step too short
    // to make any progress.
    if (this->residual > this->previousResidual * (1.0 + RSolverFluid::relaxationRiseTolerance))
    {
        this->relaxation = std::max(this->relaxation * RSolverFluid::relaxationCutFactor,
                                    RSolverFluid::minRelaxation);
    }
    else
    {
        this->relaxation = std::min(this->relaxation * RSolverFluid::relaxationGrowFactor,1.0);
    }
}

double RSolverFluid::findTimeScale() const
{
    double v = RSolverFluid::computeStreamVelocity(*this->pModel,this->nodeVelocity,true);
    double l = 1.0 / this->scales.getMetre();

    return v/l;
}

double RSolverFluid::findReScale() const
{
    double v = RSolverFluid::computeStreamVelocity(*this->pModel,this->nodeVelocity,true);
    double l = 1.0 / this->scales.getMetre();

    double Re = (this->avgU == 0.0) ? 1.0 : this->avgRo * v * l / this->avgU;

    return (Re == 0.0) ? 1.0 : 1.0e-2 / Re;
}

double RSolverFluid::findWeightScale() const
{
    double v = RSolverFluid::computeStreamVelocity(*this->pModel,this->nodeVelocity,true);
    double l = 1.0 / this->scales.getMetre();

    double ws = (this->avgU == 0.0) ? 1.0 : 1.0 *  v / (this->avgU * l * l);
    if (this->pModel->getTimeSolver().getEnabled())
    {
        ws /= this->pModel->getTimeSolver().getCurrentTimeStepSize();
    }

    return ws;
}

void RSolverFluid::computeElementScales()
{
    this->elementScales.resize(this->pModel->getNElements());
    this->elementScales.fill(0.0);

#pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->pModel->getNElements());i++)
    {
        if (!this->computableElements[uint(i)])
        {
            continue;
        }
        const RElement &rElement = this->pModel->getElement(uint(i));
        double volume = 0.0;
        if (!rElement.findVolume(this->pModel->getNodes(),volume))
        {
            continue;
        }
        if (rElement.getType() == R_ELEMENT_TETRA1)
        {
            this->elementScales[uint(i)] = std::cbrt(6.0 * volume / RConstants::pi);
        }
        else if (rElement.getType() == R_ELEMENT_HEXA1)
        {
            this->elementScales[uint(i)] = std::cbrt(volume);
        }
        else
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Failed to calculate element scales. Unsupported elemenent type \'%s\'.",
                         RElement::getName(rElement.getType()).toUtf8().constData());
        }
    }
}

void RSolverFluid::computeElementFreePressure(RRVector &values, RBVector &setValues)
{
    values.resize(this->pModel->getNElements());
    values.fill(0.0);
    setValues.resize(this->pModel->getNElements());
    setValues.fill(false);

#pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->pModel->getNSurfaces());i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(uint(i));
        for (uint j=0;j<rSurface.getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = rSurface.getBoundaryCondition(j);
            if (bc.getType() != R_BOUNDARY_CONDITION_PRESSURE_IMPLICIT)
            {
                continue;
            }
            double p = 0.0;
            uint cpos = bc.findComponentPosition(R_VARIABLE_PRESSURE);
            if (cpos != RConstants::eod)
            {
                p = bc.getComponent(cpos).get(this->pModel->getTimeSolver().getCurrentTime());
            }

            for (uint k=0;k<rSurface.size();k++)
            {
                values[rSurface.get(k)] = p;
                setValues[rSurface.get(k)] = true;
            }
            break;
        }
    }
}

void RSolverFluid::buildSparseMatrixPattern(const RBVector &elementFreePressureSetValues)
{
    this->A.clear();
    this->A.setNRows(this->nodeBook.getNEnabled());
    this->A.reserveNColumns(100);

    const uint ne = this->pModel->getNElements();
    this->elementActiveDofs.assign(ne,std::vector<uint>());
    this->elementMatrixPositions.assign(ne,std::vector<uint>());
    this->elementVectorPositions.assign(ne,std::vector<uint>());

    for (uint elementID=0;elementID<ne;elementID++)
    {
        const RElement &rElement = this->pModel->getElement(elementID);

        if (R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
        {
            if (!this->computableElements[elementID])
            {
                continue;
            }
        }
        else if (R_ELEMENT_TYPE_IS_SURFACE(rElement.getType()))
        {
            if (!elementFreePressureSetValues[elementID])
            {
                continue;
            }
        }
        else
        {
            continue;
        }

        const uint dims = 4;
        const uint ndofs = rElement.size() * dims;
        std::vector<uint> &activeDofs = this->elementActiveDofs[elementID];
        std::vector<uint> &vectorPositions = this->elementVectorPositions[elementID];
        std::vector<uint> &matrixPositions = this->elementMatrixPositions[elementID];
        activeDofs.reserve(ndofs);
        vectorPositions.assign(ndofs,RConstants::eod);
        matrixPositions.assign(ndofs*ndofs,RConstants::eod);

        for (uint m=0;m<rElement.size();m++)
        {
            uint mDims = dims*m;
            uint mIdDims = dims*rElement.getNodeId(m);
            for (uint i=0;i<dims;i++)
            {
                uint mp = 0;
                uint row = mDims+i;
                if (this->nodeBook.getValue(mIdDims+i,mp))
                {
                    vectorPositions[row] = mp;
                    activeDofs.push_back(row);
                }
            }
        }

        for (uint row : activeDofs)
        {
            uint mp = vectorPositions[row];
            for (uint column : activeDofs)
            {
                this->A.addValue(mp,vectorPositions[column],0.0);
            }
        }
    }

    for (uint elementID=0;elementID<ne;elementID++)
    {
        const std::vector<uint> &activeDofs = this->elementActiveDofs[elementID];
        std::vector<uint> &vectorPositions = this->elementVectorPositions[elementID];
        std::vector<uint> &matrixPositions = this->elementMatrixPositions[elementID];
        if (activeDofs.empty())
        {
            continue;
        }

        const uint ndofs = vectorPositions.size();
        for (uint row : activeDofs)
        {
            uint mp = vectorPositions[row];
            for (uint column : activeDofs)
            {
                uint np = vectorPositions[column];
                uint columnPosition = 0;
                R_ERROR_ASSERT(this->A.findColumnPosition(mp,np,columnPosition));
                matrixPositions[row*ndofs + column] = columnPosition;
            }
        }
    }
}

void RSolverFluid::assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    const std::vector<uint> &activeDofs = this->elementActiveDofs[elementID];
    const std::vector<uint> &vectorPositions = this->elementVectorPositions[elementID];
    const std::vector<uint> &matrixPositions = this->elementMatrixPositions[elementID];

    // Assembly final matrix system
    uint dims = 4;
    uint ndofs = rElement.size() * dims;
    for (uint row : activeDofs)
    {
        uint mp = vectorPositions[row];
        bp[mp] += fe[row];
        for (uint column : activeDofs)
        {
            uint columnPosition = matrixPositions[row*ndofs + column];
            Ap.addValueAtPosition(mp,columnPosition,Ae[row][column]);
        }
    }
}

void RSolverFluid::applyLocalRotations(unsigned int elementID, RRMatrix &Ae)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    RRMatrix T;

    bool first = true;

    for (uint i=0;i<rElement.size();i++)
    {
        uint nodeId = rElement.getNodeId(i);
        if (this->localRotations[nodeId].isActive())
        {
            if (first)
            {
                T.setIdentity(Ae.getNRows());
                first = false;
            }
            T.setBlock(this->localRotations[nodeId].getR(),4*i,4*i);
        }
    }
    if (!first)
    {
        RRMatrix Tt(T);
        Tt.transpose();

        RRMatrix Aetmp;

        RRMatrix::mlt(Tt,Ae,Aetmp);
        RRMatrix::mlt(Aetmp,T,Ae);
    }
}

double RSolverFluid::computeStreamVelocity(const RModel &rModel, const RSolverCartesianVector<RRVector> &nodeVelocity, bool averageBased)
{
    double velocity = 0.0;

    if (averageBased)
    {
#pragma omp parallel for default(shared) reduction(+:velocity)
        for (int64_t i=0;i<int64_t(rModel.getNNodes());i++)
        {
            double vx = nodeVelocity.x[uint(i)];
            double vy = nodeVelocity.y[uint(i)];
            double vz = nodeVelocity.z[uint(i)];
            velocity += std::sqrt(vx*vx + vy*vy + vz*vz);
        }
        if (rModel.getNNodes())
        {
            velocity /= double(rModel.getNNodes());
        }
    }
    else
    {
        double totalArea = 0.0;
        double totalVolurate = 0.0;

        for (uint i=0;i<rModel.getNSurfaces();i++)
        {
            const RSurface &rSurface = rModel.getSurface(i);

            bool hasVelocity = rSurface.hasBoundaryCondition(R_BOUNDARY_CONDITION_INFLOW_VELOCITY);
            bool hasVolurate = rSurface.hasBoundaryCondition(R_BOUNDARY_CONDITION_INFLOW_VOLURATE);

            double v = 0.0;
            double q = 0.0;

            if (!hasVelocity && !hasVolurate)
            {
                continue;
            }

            if (hasVelocity)
            {
                const RBoundaryCondition &rBoundaryCondition = rSurface.getBoundaryCondition(R_BOUNDARY_CONDITION_INFLOW_VELOCITY);
                uint vp = rBoundaryCondition.findComponentPosition(R_VARIABLE_VELOCITY);
                if (vp != RConstants::eod)
                {
                    v = rBoundaryCondition.getComponent(vp).get(rModel.getTimeSolver().getCurrentTime());
                }
            }

            if (hasVolurate)
            {
                const RBoundaryCondition &rBoundaryCondition = rSurface.getBoundaryCondition(R_BOUNDARY_CONDITION_INFLOW_VOLURATE);
                uint vp = rBoundaryCondition.findComponentPosition(R_VARIABLE_VOLUME_FLOW_RATE);
                if (vp != RConstants::eod)
                {
                    q = rBoundaryCondition.getComponent(vp).get(rModel.getTimeSolver().getCurrentTime());
                }
            }

            double area = rSurface.findArea(rModel.getNodes(),rModel.getElements());
            if (area < RConstants::eps)
            {
                continue;
            }
            if (hasVelocity)
            {
                q = v * area;
            }
            totalArea += area;
            totalVolurate += std::abs(q);
        }

        if (totalArea < RConstants::eps)
        {
            velocity = 1.0;
        }
        else
        {
            velocity = totalVolurate / totalArea;
        }
    }

    return (velocity < RConstants::eps) ? 1.0 : velocity;
}

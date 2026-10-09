#include <atomic>

#include <omp.h>

#include "rsolverfluid.h"
#include "rsolverfluidparticle.h"
#include "rmatrixsolver.h"

class FluidParticleMatrixContainer
{
    public:

        bool initialized;

        // Element level matricies
        RRMatrix me;    // m
        RRMatrix ce;    // c
        RRMatrix ke;    // k
        RRMatrix cte;   // c~ (SUPG weighted mass)
        RRMatrix kte;   // k~ (SUPG streamline diffusion)
        RRMatrix re;    // r (saturation limited source)
        RRMatrix vb;    // v_k . grad(N_n) with nodal velocity v_k

        // Element level vectors
        RRVector fv;    // f
        RRVector vdiv;
        RRVector sv;    // node saturation coefficients

    public:

        FluidParticleMatrixContainer() : initialized(false)
        {

        }

        void resize(uint nen)
        {
            this->me.resize(nen,nen,0.0);
            this->ce.resize(nen,nen,0.0);
            this->ke.resize(nen,nen,0.0);
            this->cte.resize(nen,nen,0.0);
            this->kte.resize(nen,nen,0.0);
            this->re.resize(nen,nen,0.0);
            this->vb.resize(nen,nen,0.0);

            this->fv.resize(nen,0.0);
            this->vdiv.resize(nen,0.0);
            this->sv.resize(nen,0.0);

            this->initialized = true;
        }

        void clear()
        {
            this->me.fill(0.0);
            this->ce.fill(0.0);
            this->ke.fill(0.0);
            this->cte.fill(0.0);
            this->kte.fill(0.0);
            this->re.fill(0.0);
            this->vb.fill(0.0);

            this->fv.fill(0.0);
            this->vdiv.fill(0.0);
            this->sv.fill(0.0);
        }
};

//! Compute SUPG stabilization parameter (Shakib, Tezduyar)
//!   tau = [ (2/dt)^2 + (2|v|/h)^2 + C (4k/h^2)^2 ]^(-1/2), C = 9 for linear elements
//! Time step term keeps tau <= dt/2. Without it tau = h/(2|v|) grows without bound as the element
//! Courant number |v|dt/h drops, while the SUPG weighted mass term tau*(v.grad(N))*N stays at full
//! strength regardless of |v|. In slow and recirculating flow this destabilizes the solution.
//! Diffusion term with C = 9 gives tau -> h^2/(12k), the limit of the nodally exact 1D parameter
//! (h/(2|v|)) * (coth(Pe) - 1/Pe), Pe = |v|h/(2k), when diffusion dominates.
//! Zero time step size means steady problem.
static double computeSupgParameter(double h, double mvh, double k, double dt)
{
    if (mvh <= 0.0 || h <= 0.0)
    {
        // No streamline direction - no streamline upwinding.
        return 0.0;
    }

    double tauInv2 = std::pow(2.0 * mvh / h, 2);
    if (dt > 0.0)
    {
        tauInv2 += std::pow(2.0 / dt, 2);
    }
    if (k > 0.0)
    {
        tauInv2 += 9.0 * std::pow(4.0 * k / (h * h), 2);
    }

    return 1.0 / std::sqrt(tauInv2);
}

RSolverFluidParticle::RSolverFluidParticle(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , streamVelocity(1.0)
    , maximumSaturation(0.0)
    , cvgC(0.0)
    , statsCounter(0)
    , statsOldResidual(0.0)
{
    this->problemType = R_PROBLEM_FLUID_PARTICLE;
}

RSolverFluidParticle::~RSolverFluidParticle()
{
    this->clearShapeDerivatives();
}

bool RSolverFluidParticle::hasConverged() const
{
    return true;
}

void RSolverFluidParticle::generateNodeRateVector()
{
    RBVector rateSetValues;
    this->generateVariableVector(R_VARIABLE_PARTICLE_RATE,this->elementRate,rateSetValues,true,this->firstRun,this->firstRun);

    this->nodeRate.fill(0.0); // Particle rate on node is meant as an input - needs to be cleared
    this->pModel->convertElementToNodeVector(this->elementRate,rateSetValues,this->nodeRate,true);
}

double RSolverFluidParticle::computeNodeSaturationCoefficient(uint nodeID) const
{
    // Only emission (positive rate) is limited - a sink must keep removing contaminant regardless of saturation.
    if (this->maximumSaturation <= 0.0 || this->nodeRate[nodeID] <= 0.0)
    {
        return 0.0;
    }
    return this->nodeRate[nodeID] / this->maximumSaturation;
}

void RSolverFluidParticle::checkMaximumSaturation() const
{
    if (this->maximumSaturation <= 0.0)
    {
        return;
    }

    uint nExceeded = 0;
    double cMax = 0.0;
    for (uint i=0;i<this->nodeConcentration.size();i++)
    {
        // Check only prescribed values - boundary conditions and on the first run also initial conditions.
        uint position = 0;
        if (!this->firstRun && this->nodeBook.getValue(i,position))
        {
            continue;
        }
        if (this->nodeConcentration[i] > this->maximumSaturation)
        {
            nExceeded++;
            cMax = std::max(cMax,this->nodeConcentration[i]);
        }
    }

    if (nExceeded > 0)
    {
        double scale = this->scales.findScaleFactor(R_VARIABLE_PARTICLE_CONCENTRATION);
        RLogger::warning("Prescribed particle concentration on %u nodes exceeds maximum saturation (%g > %g [kg/m^3]).\n",
                         nExceeded,
                         cMax/scale,
                         this->maximumSaturation/scale);
    }
}

void RSolverFluidParticle::initialize()
{
}

void RSolverFluidParticle::updateScales()
{

}

void RSolverFluidParticle::recover()
{
    this->recoveryStopWatch.reset();
    this->recoveryStopWatch.resume();
    this->recoverVariable(R_VARIABLE_PARTICLE_CONCENTRATION,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeConcentration,
                          0.0);
    this->recoverVariable(R_VARIABLE_PARTICLE_RATE,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementRate,
                          0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
    this->recoveryStopWatch.pause();
}

void RSolverFluidParticle::prepare()
{
    RLogger::info("Building matrix system\n");
    RLogger::indent();

    this->buildStopWatch.reset();
    this->assemblyStopWatch.reset();

    RBVector concentrationSetValues;

    if (this->taskIteration == 0 || this->meshChanged)
    {
        this->generateNodeBook(R_PROBLEM_FLUID_PARTICLE);
        this->findBoundaryFaces();
    }

    this->generateVariableVector(R_VARIABLE_PARTICLE_CONCENTRATION,this->elementConcentration,concentrationSetValues,true,this->firstRun,this->firstRun);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);
    // Model is downscaled at this point - diffusion coefficient [m^2/s] must be scaled as well.
    double diffusionCoefficient = std::max(this->pModel->getProblemSetup().getParticleSetup().getDiffusionCoefficient(),0.0)
                                * std::pow(this->scales.getMetre(),2) / this->scales.getSecond();
    this->elementDiffusion.resize(this->pModel->getNElements());
    this->elementDiffusion.fill(diffusionCoefficient);

    this->generateNodeRateVector();

    // Model is downscaled at this point - saturation must be scaled same as concentration.
    this->maximumSaturation = std::max(this->pModel->getProblemSetup().getParticleSetup().getMaximumSaturation(),0.0)
                            * this->scales.findScaleFactor(R_VARIABLE_PARTICLE_CONCENTRATION);

    this->pModel->convertElementToNodeVector(this->elementConcentration,concentrationSetValues,this->nodeConcentration,true);

    if (this->taskIteration == 0)
    {
        this->checkMaximumSaturation();
    }

    this->pModel->convertNodeToElementVector(this->nodeVelocity.x,this->elementVelocity.x);
    this->pModel->convertNodeToElementVector(this->nodeVelocity.y,this->elementVelocity.y);
    this->pModel->convertNodeToElementVector(this->nodeVelocity.z,this->elementVelocity.z);

    if (this->meshChanged)
    {
        this->clearShapeDerivatives();
        this->computeShapeDerivatives();
    }
    if (this->taskIteration == 0)
    {
        this->streamVelocity = RSolverFluid::computeStreamVelocity(*this->pModel,this->nodeVelocity,false);
    }

    uint nEnabled = this->nodeBook.getNEnabled();

    this->A.clear();
    this->A.setNRows(nEnabled);
    this->A.reserveNColumns(30);
    this->b.resize(nEnabled);
    this->x.resize(nEnabled);
    this->b.fill(0.0);
    this->x.fill(0.0);

    int np = omp_get_max_threads();

    QVector<RSparseMatrix> Ap;
    Ap.resize(np);
    QVector<RRVector> bp;
    bp.resize(np);

    for (int i=0;i<np;i++)
    {
        Ap[i].setNRows(nEnabled);
        Ap[i].reserveNColumns(30);
        bp[i].resize(nEnabled);
        bp[i].fill(0.0);
    }

    std::atomic<bool> abort{false};

    RMatrixManager<FluidParticleMatrixContainer> matrixManager;

    this->buildStopWatch.resume();

    // Compute element matrices
    #pragma omp parallel for default(shared) private(matrixManager)
    for (int64_t i=0;i<int64_t(this->pModel->getNElements());i++)
    {
        uint elementID = i;

        const RElement &element = this->pModel->getElement(elementID);

        if (abort.load(std::memory_order_relaxed))
        {
            continue;
        }
        try
        {
            uint nen = element.size();

            RRMatrix Ae(nen,nen,0.0);
            RRVector be(nen,0.0);

            if (R_ELEMENT_TYPE_IS_VOLUME(element.getType()))
            {
                if (!this->computableElements[elementID])
                {
                    continue;
                }
                this->computeElement(elementID,Ae,be,matrixManager);
            }
            else if (R_ELEMENT_TYPE_IS_SURFACE(element.getType()) && this->boundaryFaces[elementID])
            {
                this->computeBackflowElement(elementID,Ae,be);
            }
            this->assemblyMatrix(elementID,Ae,be,Ap[omp_get_thread_num()],bp[omp_get_thread_num()]);
        }
        catch (const RError &rError)
        {
            #pragma omp critical
            {
                RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
                abort = true;
            }
        }
    }

    this->buildStopWatch.pause();
    this->assemblyStopWatch.resume();

#pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->A.getNRows());i++)
    {
        for (int j=0;j<np;j++)
        {
            this->A.getVector(uint(i)).addVector(Ap[j].getVector(uint(i)));
            this->b[uint(i)] += bp[j][uint(i)];
        }
    }

    this->assemblyStopWatch.pause();

    if (abort)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
    }

    RLogger::unindent();
}

void RSolverFluidParticle::solve()
{
    RLogger::info("Solving matrix system\n");
    RLogger::indent();

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

    this->nodeConcentration.resize(this->pModel->getNNodes(),0.0);

    this->updateStopWatch.reset();
    this->updateStopWatch.resume();

    double cOld = RRVector::euclideanNorm(this->nodeConcentration);

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position = 0;
        if (this->nodeBook.getValue(i,position))
        {
            // Not clipped to zero: clipping undershoots adds contaminant mass which accumulates over time-steps.
            this->nodeConcentration[i] = this->x[position];
        }
    }

    double c = RRVector::euclideanNorm(this->nodeConcentration);

    this->cvgC = (c - cOld) / this->scales.findScaleFactor(R_VARIABLE_PARTICLE_CONCENTRATION);

    this->updateStopWatch.pause();

    RLogger::unindent();
}

void RSolverFluidParticle::process()
{

}

void RSolverFluidParticle::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Concentration
    uint concentrationPos = this->pModel->findVariable(R_VARIABLE_PARTICLE_CONCENTRATION);
    if (concentrationPos == RConstants::eod)
    {
        concentrationPos = this->pModel->addVariable(R_VARIABLE_PARTICLE_CONCENTRATION);

        this->pModel->getVariable(concentrationPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeConcentration),
                    RStatistics::findMaximumValue(this->nodeConcentration));
    }
    RVariable &concentration =  this->pModel->getVariable(concentrationPos);

    concentration.setApplyType(R_VARIABLE_APPLY_NODE);
    concentration.resize(1,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        concentration.setValue(0,i,this->nodeConcentration[i]);
    }

    // Relative saturation
    uint relativeSaturationPos = this->pModel->findVariable(R_VARIABLE_RELATIVE_SATURATION);
    if (this->maximumSaturation > 0.0)
    {
        if (relativeSaturationPos == RConstants::eod)
        {
            relativeSaturationPos = this->pModel->addVariable(R_VARIABLE_RELATIVE_SATURATION);

            this->pModel->getVariable(relativeSaturationPos).getVariableData().setMinMaxDisplayValue(0.0,1.0);
        }
        RVariable &relativeSaturation = this->pModel->getVariable(relativeSaturationPos);

        relativeSaturation.setApplyType(R_VARIABLE_APPLY_NODE);
        relativeSaturation.resize(1,this->pModel->getNNodes());
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            relativeSaturation.setValue(0,i,this->nodeConcentration[i] / this->maximumSaturation);
        }
    }
    else if (relativeSaturationPos != RConstants::eod)
    {
        // Saturation is no longer limited - remove stale results.
        this->pModel->removeVariable(relativeSaturationPos);
    }

    RLogger::unindent();
}

void RSolverFluidParticle::statistics()
{
    double scale = std::pow(this->scales.getSecond(),2) / this->scales.getKilogram();
    double residual = RRVector::euclideanNorm(this->b)*scale;
    double convergence = residual - this->statsOldResidual;
    this->statsOldResidual = residual;

    std::vector<RIterationInfoValue> cvgValues;
    cvgValues.push_back(RIterationInfoValue("Solver residual",residual));
    cvgValues.push_back(RIterationInfoValue("Solver convergence",convergence));
    cvgValues.push_back(RIterationInfoValue("Concentration convergence",this->cvgC));

    RIterationInfo::writeToFile(this->convergenceFileName,this->statsCounter,cvgValues);

    this->printStats(R_VARIABLE_PARTICLE_CONCENTRATION);
    this->processMonitoringPoints();

    RLogger::info("Convergence:   %-13g\n",residual);

    RLogger::info("Build time:    %9u [ms]\n",this->buildStopWatch.getMiliSeconds());
    RLogger::info("Assembly time: %9u [ms]\n",this->assemblyStopWatch.getMiliSeconds());
    RLogger::info("Solver time:   %9u [ms]\n",this->solverStopWatch.getMiliSeconds());
    RLogger::info("Update time:   %9u [ms]\n",this->updateStopWatch.getMiliSeconds());

    this->statsCounter++;
}

void RSolverFluidParticle::computeShapeDerivatives()
{
    this->shapeDerivations.resize(this->pModel->getNElements(),0);

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        uint elementID = i;

        const RElement &rElement = this->pModel->getElement(elementID);
        if (R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
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

void RSolverFluidParticle::clearShapeDerivatives()
{
    for (uint i=0;i<this->shapeDerivations.size();i++)
    {
        delete this->shapeDerivations[i];
    }
    this->shapeDerivations.clear();
}

void RSolverFluidParticle::findBoundaryFaces()
{
    this->boundaryFaces.resize(this->pModel->getNElements());
    this->boundaryFaces.fill(false);

    // Node -> computable volume elements.
    std::vector<std::vector<uint>> nodeToVolumeElements(this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        const RElement &rElement = this->pModel->getElement(i);
        if (!R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()) || !this->computableElements[i])
        {
            continue;
        }
        for (uint j=0;j<rElement.size();j++)
        {
            nodeToVolumeElements[rElement.getNodeId(j)].push_back(i);
        }
    }

    // A face of exactly one volume element lies on the boundary - faces between two volumes are skipped.
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        const RElement &rElement = this->pModel->getElement(i);
        if (!R_ELEMENT_TYPE_IS_SURFACE(rElement.getType()) || rElement.size() == 0)
        {
            continue;
        }
        uint nVolumeElements = 0;
        for (uint volumeElementID : nodeToVolumeElements[rElement.getNodeId(0)])
        {
            const RElement &rVolumeElement = this->pModel->getElement(volumeElementID);
            bool hasAllNodes = true;
            for (uint j=1;j<rElement.size() && hasAllNodes;j++)
            {
                hasAllNodes = rVolumeElement.hasNodeId(rElement.getNodeId(j));
            }
            if (hasAllNodes)
            {
                nVolumeElements++;
            }
        }
        this->boundaryFaces[i] = (nVolumeElements == 1);
    }
}

void RSolverFluidParticle::computeBackflowElement(unsigned int elementID, RRMatrix &Ae, RRVector &be)
{
    // Flow entering through a free boundary - typically backflow at a pressure outlet - carries no
    // prescribed concentration. The Galerkin advection term then feeds energy into the solution there,
    // -0.5 * int( v.n * C^2 ) > 0 for v.n < 0, and it grows without bound. Such flow is treated as an
    // inflow of clean fluid by the weak term
    //   int( max(-v.n,0) * C * w )
    // which turns that boundary contribution into dissipation. Outflow and walls (v.n >= 0) are not
    // affected. Nodes with prescribed concentration are eliminated during assembly.
    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();

    Ae.fill(0.0);
    be.fill(0.0);

    double area = 0.0;
    RR3Vector n;
    if (!element.findArea(this->pModel->getNodes(),area) ||
        !element.findNormal(this->pModel->getNodes(),n[0],n[1],n[2]))
    {
        return;
    }
    // Outward unit normal.
    n.normalize();
    if (this->inwardElements[elementID])
    {
        n[0] = -n[0];
        n[1] = -n[1];
        n[2] = -n[2];
    }

    bool unsteady = (this->pModel->getTimeSolver().getEnabled());
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    for (uint m=0;m<nen;m++)
    {
        uint nodeID = element.getNodeId(m);
        double vn = this->nodeVelocity.x[nodeID] * n[0]
                  + this->nodeVelocity.y[nodeID] * n[1]
                  + this->nodeVelocity.z[nodeID] * n[2];
        // Lumped, so that the term only damps.
        double r = std::max(-vn,0.0) * area / double(nen);
        if (unsteady)
        {
            Ae[m][m] = alpha * dt * r;
            be[m] = -(1.0 - alpha) * dt * r * this->nodeConcentration[nodeID];
        }
        else
        {
            Ae[m][m] = r;
        }
    }
}

void RSolverFluidParticle::computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager)
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

void RSolverFluidParticle::computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();
    uint nInp = RElement::getNIntegrationPoints(element.getType());

    double ro = this->elementDensity[elementID];
    double k = this->elementDiffusion[elementID];

    double ca = ro;
    ca = 1.0; // Seems like density should not be part of calculation

    Ae.fill(0.0);
    be.fill(0.0);

    FluidParticleMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());
    matrixCotainer.clear();

    // Element level matricies
    RRMatrix &me = matrixCotainer.me;
    RRMatrix &ce = matrixCotainer.ce;
    RRMatrix &ke = matrixCotainer.ke;
    RRMatrix &cte = matrixCotainer.cte;
    RRMatrix &kte = matrixCotainer.kte;
    RRMatrix &re = matrixCotainer.re;

    // Element level vectors
    RRVector &fv = matrixCotainer.fv;

    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    // Element level input -------------------------------------------
    RR3Vector ve(this->elementVelocity.x[elementID],
                 this->elementVelocity.y[elementID],
                 this->elementVelocity.z[elementID]);
    // element level velocity magnitude
    double mvh = ve.length();
    // element level velocity direction
    RR3Vector s(ve);
    s.normalize();
    // node level saturation coefficients
    RRVector &sv = matrixCotainer.sv;
    for (uint n=0;n<nen;n++)
    {
        sv[n] = this->computeNodeSaturationCoefficient(element.getNodeId(n));
    }

    for (uint intPoint=0;intPoint<nInp;intPoint++)
    {
        const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),intPoint);
        const RRVector &N = shapeFunc.getN();
        const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(intPoint);
        double detJ = this->shapeDerivations[elementID]->getJacobian(intPoint);

        double integValue = detJ * shapeFunc.getW();

        // Continuous (interpolated) velocity and its divergence at integration point for the Galerkin
        // advection term in skew-symmetric form, see computeElementConstantDerivative().
        RR3Vector vip(0.0,0.0,0.0);
        double divV = 0.0;
        for (uint n=0;n<nen;n++)
        {
            uint nodeID = element.getNodeId(n);
            RR3Vector vn(this->nodeVelocity.x[nodeID],this->nodeVelocity.y[nodeID],this->nodeVelocity.z[nodeID]);
            vip[0] += N[n] * vn[0];
            vip[1] += N[n] * vn[1];
            vip[2] += N[n] * vn[2];
            divV += B[n][0] * vn[0] + B[n][1] * vn[1] + B[n][2] * vn[2];
        }

        // advection operator (v . grad N)
        RRVector &vdiv = matrixCotainer.vdiv;
        vdiv.fill(0.0);
        // element length scale
        double h = 0.0;

        for (uint m=0;m<nen;m++)
        {
            vdiv[m] += ve[0] * B[m][0] + ve[1] * B[m][1] + ve[2] * B[m][2];
            h += std::fabs(s[0]*B[m][0] + s[1]*B[m][1] + s[2]*B[m][2]);
        }
        if (h != 0.0)
        {
            h = 2.0/h;
        }
        // SUPG stabilization parameter
        double Tsupg = computeSupgParameter(h,mvh,k,unsteady ? dt : 0.0);

        // Particle rate at integration point
        double fip = 0.0;
        // Saturation coefficient at integration point
        double sip = 0.0;
        for (uint n=0;n<nen;n++)
        {
            fip += this->nodeRate[element.getNodeId(n)] * N[n];
            sip += sv[n] * N[n];
        }

        for (uint m=0;m<nen;m++)
        {
            for (uint n=0;n<nen;n++)
            {
                // m matrix
                me[m][n] = unsteady ? ca * N[m] * N[n] : 0.0;
                // c matrix (skew-symmetric form)
                ce[m][n] = ca * N[m] * (  vip[0] * B[n][0] + vip[1] * B[n][1] + vip[2] * B[n][2]
                                        + 0.5 * divV * N[n]);
                // k matrix
                ke[m][n] = k * (B[m][0] * B[n][0] + B[m][1] * B[n][1] + B[m][2] * B[n][2]);
                // c~ matrix
                cte[m][n] = unsteady ? Tsupg * ca * vdiv[m] * N[n] : 0.0;
                // k~ matrix
                kte[m][n] = Tsupg * ca * vdiv[m] * vdiv[n];
                // r matrix (SUPG weighted saturation limit, Galerkin part is added below)
                re[m][n] = Tsupg * vdiv[m] * sip * N[n];
            }
            // r matrix (Galerkin saturation limit lumped same as the source so that it cancels the source at saturation)
            re[m][m] += sv[m] * N[m];
            // f vector (Galerkin + SUPG weighted source)
            fv[m] = this->nodeRate[element.getNodeId(m)] * N[m] + Tsupg * vdiv[m] * fip;
        }

        // Assembly element level matrixes
        for (uint m=0;m<nen;m++)
        {
            for (uint n=0;n<nen;n++)
            {
                if (unsteady)
                {
                    Ae[m][n] += integValue * (  me[m][n] + cte[m][n]
                                              + alpha * dt * (ce[m][n] + ke[m][n] + kte[m][n] + re[m][n]));
                }
                else
                {
                    Ae[m][n] += integValue * (ce[m][n] + ke[m][n] + kte[m][n] + re[m][n]);
                }
            }
            if (unsteady)
            {
                be[m] += integValue * dt * fv[m];
                for (uint n=0;n<nen;n++)
                {
                    be[m] += integValue
                           * (me[m][n] + cte[m][n] - (1.0 - alpha) * dt * (ce[m][n] + ke[m][n] + kte[m][n] + re[m][n]))
                           * this->nodeConcentration[element.getNodeId(n)];
                }
            }
            else
            {
                be[m] += integValue * fv[m];
            }
        }
    }
}

void RSolverFluidParticle::computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidParticleMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();

    double ro = this->elementDensity[elementID];
    double k = this->elementDiffusion[elementID];

    double ca = ro;
    ca = 1.0; // Seems like density should not be part of calculation

    Ae.fill(0.0);
    be.fill(0.0);

    FluidParticleMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());
    matrixCotainer.clear();

    // Element level matricies
    RRMatrix &me = matrixCotainer.me;
    RRMatrix &ce = matrixCotainer.ce;
    RRMatrix &ke = matrixCotainer.ke;
    RRMatrix &cte = matrixCotainer.cte;
    RRMatrix &kte = matrixCotainer.kte;
    RRMatrix &re = matrixCotainer.re;

    // Element level vectors
    RRVector &fv = matrixCotainer.fv;

    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    // Element level input -------------------------------------------
    RR3Vector ve(this->elementVelocity.x[elementID],
                 this->elementVelocity.y[elementID],
                 this->elementVelocity.z[elementID]);
    // element level velocity magnitude
    double mvh = ve.length();
    // element level velocity direction
    RR3Vector s(ve);
    s.normalize();
    // node level saturation coefficients
    RRVector &sv = matrixCotainer.sv;
    for (uint n=0;n<nen;n++)
    {
        sv[n] = this->computeNodeSaturationCoefficient(element.getNodeId(n));
    }

    const RRVector &iN = RElement::getMassVector(element.getType());
    const RRMatrix &iNiN = RElement::getMassMatrix(element.getType());
    double wt = RElement::getTotalWeightFactor(element.getType());
    const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(0);

    // advection operator (v . grad N)
    RRVector &vdiv = matrixCotainer.vdiv;
    vdiv.fill(0.0);
    // element length scale
    double h = 0.0;

    for (uint m=0;m<nen;m++)
    {
        vdiv[m] += ve[0] * B[m][0] + ve[1] * B[m][1] + ve[2] * B[m][2];
        h += std::fabs(s[0]*B[m][0] + s[1]*B[m][1] + s[2]*B[m][2]);
    }
    if (h != 0.0)
    {
        h = 2.0/h;
    }
    // SUPG stabilization parameter
    double Tsupg = computeSupgParameter(h,mvh,k,unsteady ? dt : 0.0);

    // Galerkin advection term in skew-symmetric form with continuous (nodal) velocity,
    //   int( N_m * ( v . grad(C) + 0.5 * div(v) * C ) )
    // Discrete velocity of the flow solution is not exactly divergence free and advection with an
    // element-constant velocity, whose normal component jumps across element faces, can create
    // energy - with little SUPG diffusion (small time step) the solution then grows without bound.
    // In this form the face contributions of a continuous velocity cancel, so advection only moves
    // the contaminant. For a divergence free velocity the added term vanishes.
    RRMatrix &vb = matrixCotainer.vb;
    double divV = 0.0;
    for (uint o=0;o<nen;o++)
    {
        uint nodeID = element.getNodeId(o);
        RR3Vector vo(this->nodeVelocity.x[nodeID],this->nodeVelocity.y[nodeID],this->nodeVelocity.z[nodeID]);
        for (uint n=0;n<nen;n++)
        {
            vb[o][n] = vo[0] * B[n][0] + vo[1] * B[n][1] + vo[2] * B[n][2];
        }
        divV += vb[o][o];
    }

    // Integral of particle rate over the element (reference coordinates)
    double fInt = 0.0;
    for (uint n=0;n<nen;n++)
    {
        fInt += this->nodeRate[element.getNodeId(n)] * iN[n];
    }

    for (uint m=0;m<nen;m++)
    {
        for (uint n=0;n<nen;n++)
        {
            // m matrix
            me[m][n] = unsteady ? ca * iNiN[m][n] : 0.0;
            // c matrix (skew-symmetric form)
            ce[m][n] = 0.5 * divV * iNiN[m][n];
            for (uint o=0;o<nen;o++)
            {
                ce[m][n] += iNiN[m][o] * vb[o][n];
            }
            ce[m][n] *= ca;
            // k matrix
            ke[m][n] = k * wt * (B[m][0] * B[n][0] + B[m][1] * B[n][1] + B[m][2] * B[n][2]);
            // c~ matrix
            cte[m][n] = unsteady ? Tsupg * ca * vdiv[m] * iN[n] : 0.0;
            // k~ matrix
            kte[m][n] = Tsupg * ca * vdiv[m] * vdiv[n] * wt;
            // r matrix (SUPG weighted saturation limit, Galerkin part is added below)
            re[m][n] = 0.0;
            for (uint o=0;o<nen;o++)
            {
                re[m][n] += Tsupg * vdiv[m] * sv[o] * iNiN[o][n];
            }
        }
        // r matrix (Galerkin saturation limit lumped same as the source so that it cancels the source at saturation)
        re[m][m] += sv[m] * iN[m];
        // f vector (Galerkin + SUPG weighted source)
        fv[m] = this->nodeRate[element.getNodeId(m)] * iN[m] + Tsupg * vdiv[m] * fInt;
    }

    // Assembly element level matrixes
    for (uint m=0;m<nen;m++)
    {
        for (uint n=0;n<nen;n++)
        {
            if (unsteady)
            {
                Ae[m][n] = me[m][n] + cte[m][n]
                         + alpha * dt * (ce[m][n] + ke[m][n] + kte[m][n] + re[m][n]);
            }
            else
            {
                Ae[m][n] = ce[m][n] + ke[m][n] + kte[m][n] + re[m][n];
            }
        }
        if (unsteady)
        {
            be[m] = dt * fv[m];
            for (uint n=0;n<nen;n++)
            {
                be[m] += (me[m][n] + cte[m][n] - (1.0 - alpha) * dt * (ce[m][n] + ke[m][n] + kte[m][n] + re[m][n]))
                       * this->nodeConcentration[element.getNodeId(n)];
            }
        }
        else
        {
            be[m] = fv[m];
        }
    }

    double detJ = this->shapeDerivations[elementID]->getJacobian(0);
    Ae *= detJ;
    be *= detJ;
}

void RSolverFluidParticle::assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    RRVector fe(be);

    // Apply explicit boundary conditions.
    for (uint m=0;m<rElement.size();m++)
    {
        uint mp = 0;
        if (!this->nodeBook.getValue(rElement.getNodeId(m),mp))
        {
            for (uint n=0;n<rElement.size();n++)
            {
                fe[n] -= Ae[n][m] * this->nodeConcentration[rElement.getNodeId(m)];
            }
        }
    }

    // Assembly final matrix system
    for (uint m=0;m<rElement.size();m++)
    {
        uint mp = 0;

        if (this->nodeBook.getValue(rElement.getNodeId(m),mp))
        {
            this->b[mp] += fe[m];

            for (uint n=0;n<rElement.size();n++)
            {
                uint np = 0;

                if (this->nodeBook.getValue(rElement.getNodeId(n),np))
                {
                    this->A.addValue(mp,np,Ae[m][n]);
                }
            }
        }
    }
}

void RSolverFluidParticle::assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be, RSparseMatrix &Ap, RRVector &bp)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    RRVector fe(be);

    // Apply explicit boundary conditions.
    for (uint m=0;m<rElement.size();m++)
    {
        uint mp = 0;
        uint nodeID = rElement.getNodeId(m);
        if (!this->nodeBook.getValue(nodeID,mp))
        {
            for (uint n=0;n<rElement.size();n++)
            {
                fe[n] -= Ae[n][m] * this->nodeConcentration[nodeID];
            }
        }
    }

    // Assembly final matrix system
    for (uint m=0;m<rElement.size();m++)
    {
        uint mp = 0;
        uint nodeIDm = rElement.getNodeId(m);

        if (this->nodeBook.getValue(nodeIDm,mp))
        {
            bp[mp] += fe[m];

            for (uint n=0;n<rElement.size();n++)
            {
                uint np = 0;
                uint nodeIDn = rElement.getNodeId(n);

                if (this->nodeBook.getValue(nodeIDn,np))
                {
                    Ap.addValue(mp,np,Ae[m][n]);
                }
            }
        }
    }
}

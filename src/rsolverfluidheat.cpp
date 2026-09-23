#include <algorithm>
#include <atomic>
#include <cmath>

#include <omp.h>

#include "rsolverfluid.h"
#include "rsolverfluidheat.h"
#include "rsolverheat.h"
#include "rmatrixsolver.h"

class FluidHeatMatrixContainer
{
    public:

        bool initialized;

        // Element level matricies
        RRMatrix me;    // m
        RRMatrix ce;    // c
        RRMatrix ke;    // k
        RRMatrix cte;   // c~
        RRMatrix kte;   // k~
        RRMatrix yte;   // y~

        // Element level vectors
        RRVector fv;    // f
        RRVector vdiv;

    public:

        FluidHeatMatrixContainer() : initialized(false)
        {

        }

        void resize(uint nen)
        {
            this->me.resize(nen,nen,0.0);
            this->ce.resize(nen,nen,0.0);
            this->ke.resize(nen,nen,0.0);
            this->cte.resize(nen,nen,0.0);
            this->kte.resize(nen,nen,0.0);
            this->yte.resize(nen,nen,0.0);

            this->fv.resize(nen,0.0);
            this->vdiv.resize(nen,0.0);

            this->initialized = true;
        }

        void clear()
        {
            this->me.fill(0.0);
            this->ce.fill(0.0);
            this->ke.fill(0.0);
            this->cte.fill(0.0);
            this->kte.fill(0.0);
            this->yte.fill(0.0);

            this->fv.fill(0.0);
            this->vdiv.fill(0.0);
        }
};

const QString RSolverFluidHeat::wallHeatTransferCoefficientKey("fluid-wall-heat-transfer-coefficient");
const QString RSolverFluidHeat::wallFluidTemperatureKey("fluid-wall-fluid-temperature");

RSolverFluidHeat::RSolverFluidHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , streamVelocity(1.0)
    , cvgT(0.0)
    , wallCoupled(false)
    , wallRelaxation(1.0)
    , statsCounter(0)
    , statsOldResidual(0.0)
{
    this->problemType = R_PROBLEM_FLUID_HEAT;
}

RSolverFluidHeat::~RSolverFluidHeat()
{
    this->clearShapeDerivatives();
}

bool RSolverFluidHeat::hasConverged() const
{
    // Without walls held at the solid temperature there is nothing to iterate
    // on - one solve is the answer.
    if (!this->wallCoupled)
    {
        return true;
    }
    // A task group with no convergence value runs all of its iterations.
    if (this->taskCvgValue <= 0.0)
    {
        return false;
    }
    // The first coupled pass has no previous one to compare against.
    if (this->taskIteration < 1)
    {
        return false;
    }
    return (this->cvgT < this->taskCvgValue);
}

double RSolverFluidHeat::findTemperatureScale() const
{
    return 1.0;
}

void RSolverFluidHeat::generateNodeHeatVector()
{
    RBVector heatSetValues;
    // The Heat boundary condition prescribes the total heat input for the whole
    // entity - it is spread over the entity measure to give the source density.
    this->generateHeatVector(this->elementHeat,heatSetValues);

    this->nodeHeat.fill(0.0); // Heat on node is meant as an input - needs to be cleared
    this->pModel->convertElementToNodeVector(this->elementHeat,heatSetValues,this->nodeHeat,true);

    RRVector qv(this->pModel->getNNodes(),0.0);
    RUVector qc(this->pModel->getNNodes(),0);

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        double q = this->elementJouleHeat[i] + this->elementRadiativeHeat[i];

        if (q == 0.0)
        {
            continue;
        }

        const RElement &rElement(this->pModel->getElement(i));
        for (uint j=0;j<rElement.size();j++)
        {
            qv[rElement.getNodeId(j)] += q;
            qc[rElement.getNodeId(j)]++;
        }
    }

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        if (qc[i] == 0)
        {
            continue;
        }
        this->nodeHeat[i] += qv[i] / double(qc[i]);
    }
}

void RSolverFluidHeat::findWallElements()
{
    this->wallFluidElements.resize(this->pModel->getNElements());
    this->wallFluidElements.fill(RConstants::eod);
    this->wallNodes.resize(this->pModel->getNNodes());
    this->wallNodes.fill(false);

    // The relaxation starts over with every task run - every time step.
    this->wallTemperature.clear();
    this->wallResidual.clear();
    this->wallRelaxation = 1.0;

    // Index fluid volume elements by node so the search below stays local.
    std::vector<std::vector<uint>> nodeToFluidElements(this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        const RElement &rElement = this->pModel->getElement(i);
        if (!this->computableElements[i] || !R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
        {
            continue;
        }
        for (uint j=0;j<rElement.size();j++)
        {
            nodeToFluidElements[rElement.getNodeId(j)].push_back(i);
        }
    }

    // A wall is a surface element carrying the Forced convection condition whose
    // nodes are all shared with a fluid volume element - the mesh is conformal,
    // so the fluid element behind the wall contains the whole face.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);
        if (!rSurface.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED))
        {
            continue;
        }

        for (uint j=0;j<rSurface.size();j++)
        {
            uint elementID = rSurface.get(j);
            const RElement &rElement = this->pModel->getElement(elementID);

            for (uint fluidElementID : nodeToFluidElements[rElement.getNodeId(0)])
            {
                const RElement &rFluidElement = this->pModel->getElement(fluidElementID);
                uint nNodesFound = 0;
                for (uint k=0;k<rElement.size();k++)
                {
                    if (rFluidElement.hasNodeId(rElement.getNodeId(k)))
                    {
                        nNodesFound++;
                    }
                }
                if (nNodesFound == rElement.size())
                {
                    this->wallFluidElements[elementID] = fluidElementID;
                    for (uint k=0;k<rElement.size();k++)
                    {
                        this->wallNodes[rElement.getNodeId(k)] = true;
                    }
                    break;
                }
            }
        }
    }
}

void RSolverFluidHeat::applyWallTemperature()
{
    this->wallCoupled = false;
    this->coupledWallNodes.resize(this->pModel->getNNodes());
    this->coupledWallNodes.fill(false);

    if (this->solidNodeTemperature.size() != this->pModel->getNNodes())
    {
        return;
    }

    // A node prescribed explicitly keeps its value - the solid temperature only
    // replaces the adiabatic wall.
    RBVector disabledPositions(this->nodeBook.size(),false);
    RBVector &coupledNodes = this->coupledWallNodes;
    for (uint i=0;i<this->nodeBook.size();i++)
    {
        uint position = 0;
        if (!this->nodeBook.getValue(i,position))
        {
            disabledPositions[i] = true;
        }
        else if (this->wallNodes[i])
        {
            disabledPositions[i] = true;
            coupledNodes[i] = true;
            this->wallCoupled = true;
        }
    }

    if (!this->wallCoupled)
    {
        return;
    }

    RSolverGeneric::rebuildNodeBook(this->nodeBook,disabledPositions);

    if (this->wallTemperature.size() != this->pModel->getNNodes())
    {
        // First coupled pass - take the solid temperature as it is.
        this->wallTemperature = this->solidNodeTemperature;
    }
    else
    {
        // Aitken relaxation. With the residual r = T_solid - T_wall of this pass
        // and of the previous one, the factor is updated as
        //
        //   w = -w * (r_old . (r - r_old)) / |r - r_old|^2
        //
        // which is the secant step of the fixed point iteration - exact for a
        // linear problem, where it lands on the coupled solution at once.
        RRVector residual(this->pModel->getNNodes(),0.0);
        for (uint i=0;i<residual.size();i++)
        {
            if (coupledNodes[i])
            {
                residual[i] = this->solidNodeTemperature[i] - this->wallTemperature[i];
            }
        }

        if (this->wallResidual.size() == residual.size())
        {
            double numerator = 0.0;
            double denominator = 0.0;
            for (uint i=0;i<residual.size();i++)
            {
                double dr = residual[i] - this->wallResidual[i];
                numerator += this->wallResidual[i] * dr;
                denominator += dr * dr;
            }
            if (denominator > RConstants::eps * RConstants::eps)
            {
                // Bounded, so one noisy pass can not throw the walls far off.
                const double maxRelaxation = 100.0;
                this->wallRelaxation = std::clamp(-this->wallRelaxation * numerator / denominator,-maxRelaxation,maxRelaxation);
            }
        }
        this->wallResidual = residual;

        for (uint i=0;i<residual.size();i++)
        {
            if (coupledNodes[i])
            {
                this->wallTemperature[i] = std::max(this->wallTemperature[i] + this->wallRelaxation * residual[i],0.0);
            }
        }
    }

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        if (coupledNodes[i])
        {
            this->nodeTemperature[i] = this->wallTemperature[i];
        }
    }
}

void RSolverFluidHeat::computeWallHeatTransfer()
{
    this->elementWallHtc.resize(this->pModel->getNElements());
    this->elementWallHtc.fill(-1.0);
    this->elementWallHtt.resize(this->pModel->getNElements());
    this->elementWallHtt.fill(0.0);

    // With the gradient of the fluid element g_i = dN_i/dn along the wall normal
    // pointing into the fluid, and the wall nodes all at T_w, the heat flux
    // leaving the solid is
    //
    //   q = -k * sum_i(g_i * T_i) = k * G * (T_w - T_ref)
    //
    // where the sum runs over the nodes off the wall, G = sum_i(g_i) and
    // T_ref = sum_i(g_i * T_i) / G - the derivatives of all nodes sum to zero.
    // For a linear tetrahedron G is the reciprocal of its height over the wall
    // and T_ref the temperature of the node opposite, so h = k * G is exactly the
    // conductance of the first fluid element.
    //
    // That gradient is only first order accurate in a thin boundary layer. Once
    // the walls are held at the solid temperature, the flux is instead taken
    // from the residual of the fluid system at the wall nodes - the flux a
    // single solve of both domains would see - and T_ref is set so that the same
    // conductance reproduces it:
    //
    //   T_ref = T_w - q / h
    //
    // The conductance is kept as the coefficient either way. It is what the
    // heat solver couples the wall through, and a value close to the true
    // sensitivity of the fluid flux keeps the alternation of the two solves
    // short.
    RRVector nodeFlux;
    if (this->wallCoupled)
    {
        nodeFlux = this->computeWallReaction();

        // Turn the nodal heat into a flux density with the wall area lumped to
        // the nodes.
        RRVector nodeArea(this->pModel->getNNodes(),0.0);
        for (uint i=0;i<this->pModel->getNElements();i++)
        {
            if (this->wallFluidElements[i] == RConstants::eod)
            {
                continue;
            }
            const RElement &rWallElement = this->pModel->getElement(i);
            double area = 0.0;
            if (!rWallElement.findArea(this->pModel->getNodes(),area))
            {
                continue;
            }
            for (uint j=0;j<rWallElement.size();j++)
            {
                nodeArea[rWallElement.getNodeId(j)] += area / double(rWallElement.size());
            }
        }
        for (uint i=0;i<nodeFlux.size();i++)
        {
            nodeFlux[i] = (nodeArea[i] > 0.0) ? nodeFlux[i] / nodeArea[i] : 0.0;
        }
    }

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        uint fluidElementID = this->wallFluidElements[i];
        if (fluidElementID == RConstants::eod || !this->shapeDerivations[fluidElementID])
        {
            continue;
        }

        const RElement &rWallElement = this->pModel->getElement(i);
        const RElement &rFluidElement = this->pModel->getElement(fluidElementID);

        RR3Vector normal;
        if (!rWallElement.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]))
        {
            continue;
        }

        RR3Vector wallCenter;
        RR3Vector fluidCenter;
        rWallElement.findCenter(this->pModel->getNodes(),wallCenter[0],wallCenter[1],wallCenter[2]);
        rFluidElement.findCenter(this->pModel->getNodes(),fluidCenter[0],fluidCenter[1],fluidCenter[2]);
        RR3Vector inward;
        RR3Vector::subtract(fluidCenter,wallCenter,inward);
        if (RR3Vector::dot(normal,inward) < 0.0)
        {
            normal *= -1.0;
        }

        // Derivatives averaged over the element - exact for the constant
        // derivative elements, the element centre value for the others.
        uint nInp = RElement::hasConstantDerivative(rFluidElement.getType()) ? 1 : RElement::getNIntegrationPoints(rFluidElement.getType());

        double G = 0.0;
        double GT = 0.0;
        for (uint m=0;m<rFluidElement.size();m++)
        {
            uint nodeID = rFluidElement.getNodeId(m);
            if (rWallElement.hasNodeId(nodeID))
            {
                continue;
            }
            double g = 0.0;
            for (uint k=0;k<nInp;k++)
            {
                const RRMatrix &B = this->shapeDerivations[fluidElementID]->getDerivative(k);
                g += (B[m][0]*normal[0] + B[m][1]*normal[1] + B[m][2]*normal[2]) / double(nInp);
            }
            G += g;
            GT += g * this->nodeTemperature[nodeID];
        }

        if (G <= RConstants::eps)
        {
            continue;
        }

        this->elementWallHtc[i] = this->elementConduction[fluidElementID] * G;
        this->elementWallHtt[i] = GT / G;

        if (!this->wallCoupled)
        {
            continue;
        }

        bool coupled = true;
        double wallTemperature = 0.0;
        double wallFlux = 0.0;
        for (uint j=0;j<rWallElement.size();j++)
        {
            uint nodeID = rWallElement.getNodeId(j);
            coupled = coupled && this->coupledWallNodes[nodeID];
            wallTemperature += this->nodeTemperature[nodeID] / double(rWallElement.size());
            wallFlux += nodeFlux[nodeID] / double(rWallElement.size());
        }
        if (coupled)
        {
            this->elementWallHtt[i] = wallTemperature - wallFlux / this->elementWallHtc[i];
        }
    }
}

RRVector RSolverFluidHeat::computeWallReaction()
{
    RRVector reaction(this->pModel->getNNodes(),0.0);

    bool unsteady = this->pModel->getTimeSolver().getEnabled();

    // The element matrices of a transient solve are built against the previous
    // time level, so it is put back for the rebuild.
    RRVector nodeTemperatureNew(this->nodeTemperature);
    if (unsteady && this->nodeTemperatureOld.size() == this->nodeTemperature.size())
    {
        this->nodeTemperature = this->nodeTemperatureOld;
    }

    RMatrixManager<FluidHeatMatrixContainer> matrixManager;

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        const RElement &rElement = this->pModel->getElement(i);
        if (!this->computableElements[i] || !R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
        {
            continue;
        }

        bool touchesWall = false;
        for (uint j=0;j<rElement.size() && !touchesWall;j++)
        {
            touchesWall = this->coupledWallNodes[rElement.getNodeId(j)];
        }
        if (!touchesWall)
        {
            continue;
        }

        uint nen = rElement.size();
        RRMatrix Ae(nen,nen,0.0);
        RRVector be(nen,0.0);
        this->computeElement(i,Ae,be,matrixManager);

        for (uint m=0;m<nen;m++)
        {
            uint nodeID = rElement.getNodeId(m);
            if (!this->coupledWallNodes[nodeID])
            {
                continue;
            }
            double r = -be[m];
            for (uint n=0;n<nen;n++)
            {
                r += Ae[m][n] * nodeTemperatureNew[rElement.getNodeId(n)];
            }
            reaction[nodeID] += r;
        }
    }

    this->nodeTemperature = nodeTemperatureNew;

    // A transient system is the heat balance over one time step.
    if (unsteady)
    {
        double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
        if (dt > 0.0)
        {
            reaction *= 1.0 / dt;
        }
    }

    return reaction;
}

void RSolverFluidHeat::storeSharedData()
{
    this->RSolverGeneric::storeSharedData();

    // Shared in SI units - the heat solver works in scales of its own.
    double htcScale = this->scales.findScaleFactor(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
    double temperatureScale = this->scales.findScaleFactor(R_VARIABLE_TEMPERATURE);

    RRVector htc(this->elementWallHtc);
    RRVector htt(this->elementWallHtt);
    for (uint i=0;i<htc.size();i++)
    {
        if (htc[i] >= 0.0)
        {
            htc[i] /= htcScale;
            htt[i] /= temperatureScale;
        }
    }

    this->pSharedData->addData(RSolverFluidHeat::wallHeatTransferCoefficientKey,htc);
    this->pSharedData->addData(RSolverFluidHeat::wallFluidTemperatureKey,htt);
}

void RSolverFluidHeat::recoverSharedData()
{
    this->RSolverGeneric::recoverSharedData();

    this->solidNodeTemperature.clear();
    if (this->pSharedData->hasData(RSolverHeat::solidNodeTemperatureKey,this->pModel->getNNodes()))
    {
        this->solidNodeTemperature = this->pSharedData->findData(RSolverHeat::solidNodeTemperatureKey);
        this->solidNodeTemperature *= this->scales.findScaleFactor(R_VARIABLE_TEMPERATURE);
    }
}

void RSolverFluidHeat::initialize()
{
}

void RSolverFluidHeat::updateScales()
{
    this->scales.setMetre(this->findMeshScale());
    this->scales.setKelvin(this->findTemperatureScale());
}

void RSolverFluidHeat::recover()
{
    this->recoveryStopWatch.reset();
    this->recoveryStopWatch.resume();
    this->recoverVariable(R_VARIABLE_TEMPERATURE,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeTemperature,
                          RVariable::getInitValue(R_VARIABLE_TEMPERATURE));
    this->recoverVariable(R_VARIABLE_HEAT,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementHeat,
                          0.0);
    this->recoverVariable(R_VARIABLE_HEAT_RADIATION,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementRadiativeHeat,
                          0.0);
    this->recoverVariable(R_VARIABLE_JOULE_HEAT,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementJouleHeat,
                          0.0);

    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
    this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
    this->recoveryStopWatch.pause();
}

void RSolverFluidHeat::prepare()
{
    RLogger::info("Building matrix system\n");
    RLogger::indent();

    this->buildStopWatch.reset();
    this->assemblyStopWatch.reset();

    RBVector temperatureSetValues;

    if (this->taskIteration == 0 || this->meshChanged || this->wallFluidElements.size() != this->pModel->getNElements())
    {
        this->findWallElements();
    }

    // Rebuilt on every pass - the wall nodes leave the system once the heat
    // solver has provided a solid temperature for them.
    this->generateNodeBook(R_PROBLEM_FLUID_HEAT);

    this->generateVariableVector(R_VARIABLE_TEMPERATURE,this->elementTemperature,temperatureSetValues,true,this->firstRun,this->firstRun);
    this->generateMaterialVecor(RMaterialProperty::ThermalConductivity,this->elementConduction);
    this->generateMaterialVecor(RMaterialProperty::HeatCapacity,this->elementCapacity);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);

    this->generateNodeHeatVector();

    this->pModel->convertElementToNodeVector(this->elementTemperature,temperatureSetValues,this->nodeTemperature,true);

    this->applyWallTemperature();

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

    RMatrixManager<FluidHeatMatrixContainer> matrixManager;

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

void RSolverFluidHeat::solve()
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

    this->nodeTemperature.resize(this->pModel->getNNodes(),0.0);

    this->updateStopWatch.reset();
    this->updateStopWatch.resume();

    // Kept for the wall reaction, whose element matrices are built against it.
    this->nodeTemperatureOld = this->nodeTemperature;

    // Relative size of the change - ||dT|| / ||T|| - so it can be compared
    // against the convergence value of the task group.
    double dtNorm = 0.0;
    double tNorm = 0.0;

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position = 0;
        if (this->nodeBook.getValue(i,position))
        {
            double t = std::max(this->x[position],0.0);
            double dt = t - this->nodeTemperature[i];
            dtNorm += dt*dt;
            tNorm += t*t;
            this->nodeTemperature[i] = t;
        }
    }

    this->cvgT = (tNorm > 0.0) ? std::sqrt(dtNorm/tNorm) : 0.0;

    this->updateStopWatch.pause();

    RLogger::unindent();
}

void RSolverFluidHeat::process()
{
    double Qx;
    double Qy;
    double Qz;

    // Initialize heat flux vector vector
    this->elementHeatFlux.resize(this->pModel->getNElements(),RR3Vector(0.0,0.0,0.0));

    // Process volume elements.
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        uint elementID = i;
        const RElement &element = this->pModel->getElement(elementID);
        if (!R_ELEMENT_TYPE_IS_VOLUME(element.getType()) || !this->computableElements[elementID])
        {
            continue;
        }

        uint nInp = RElement::getNIntegrationPoints(element.getType());
        RRMatrix B(element.size(),3);

        Qx = Qy = Qz = 0.0;

        B.fill(0.0);

        // Conduction
        for (uint k=0;k<nInp;k++)
        {
            const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
            const RRMatrix &dN = shapeFunc.getDN();
            RRMatrix J, Rt;
            double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

            for (uint m=0;m<dN.getNRows();m++)
            {
                B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]) * detJ / double(nInp);
                B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]) * detJ / double(nInp);
                B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]) * detJ / double(nInp);
            }
        }

        for (uint m=0;m<element.size();m++)
        {
            uint nodeID = element.getNodeId(m);

            Qx -= B[m][0] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
            Qy -= B[m][1] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
            Qz -= B[m][2] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
        }

        this->elementHeatFlux[elementID][0] = Qx;
        this->elementHeatFlux[elementID][1] = Qy;
        this->elementHeatFlux[elementID][2] = Qz;
    }

    this->computeWallHeatTransfer();
}

void RSolverFluidHeat::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Temperature
    uint temperaturePos = this->pModel->findVariable(R_VARIABLE_TEMPERATURE);
    if (temperaturePos == RConstants::eod)
    {
        temperaturePos = this->pModel->addVariable(R_VARIABLE_TEMPERATURE);
        this->pModel->getVariable(temperaturePos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeTemperature),
                    RStatistics::findMaximumValue(this->nodeTemperature));
    }
    RVariable &temperature =  this->pModel->getVariable(temperaturePos);

    temperature.setApplyType(R_VARIABLE_APPLY_NODE);
    temperature.resize(1,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        temperature.setValue(0,i,this->nodeTemperature[i]);
    }

    // Heat flux
    uint heatFluxPos = this->pModel->findVariable(R_VARIABLE_HEAT_FLUX);
    if (heatFluxPos == RConstants::eod)
    {
        heatFluxPos = this->pModel->addVariable(R_VARIABLE_HEAT_FLUX);
        this->pModel->getVariable(heatFluxPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumMagnitude(this->elementHeatFlux),
                    RStatistics::findMaximumMagnitude(this->elementHeatFlux));
    }
    RVariable &heatFlux =  this->pModel->getVariable(heatFluxPos);

    heatFlux.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    heatFlux.resize(3,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        heatFlux.setValue(0,i,this->elementHeatFlux[i][0]);
        heatFlux.setValue(1,i,this->elementHeatFlux[i][1]);
        heatFlux.setValue(2,i,this->elementHeatFlux[i][2]);
    }

    RLogger::unindent();
}

void RSolverFluidHeat::statistics()
{
    double scale = std::pow(this->scales.getSecond(),2) / this->scales.getKilogram();
    double residual = RRVector::euclideanNorm(this->b)*scale;
    double convergence = residual - this->statsOldResidual;
    this->statsOldResidual = residual;

    std::vector<RIterationInfoValue> cvgValues;
    cvgValues.push_back(RIterationInfoValue("Solver residual",residual));
    cvgValues.push_back(RIterationInfoValue("Solver convergence",convergence));
    cvgValues.push_back(RIterationInfoValue("Temperature convergence",this->cvgT));

    RIterationInfo::writeToFile(this->convergenceFileName,this->statsCounter,cvgValues);

    this->printStats(R_VARIABLE_TEMPERATURE);
    this->printStats(R_VARIABLE_HEAT_FLUX);
    this->processMonitoringPoints();

    RLogger::info("Convergence:   %-13g\n",residual);
    if (this->wallCoupled)
    {
        RLogger::info("Walls held at the solid temperature - convergence-T: %-13g\n",this->cvgT);
        if (this->taskCvgValue > 0.0)
        {
            RLogger::info("Convergence target: %-13g%s\n",
                          this->taskCvgValue,
                          this->hasConverged() ? " (reached)" : "");
        }
    }
    RLogger::info("Build time:    %9u [ms]\n",this->buildStopWatch.getMiliSeconds());
    RLogger::info("Assembly time: %9u [ms]\n",this->assemblyStopWatch.getMiliSeconds());
    RLogger::info("Solver time:   %9u [ms]\n",this->solverStopWatch.getMiliSeconds());
    RLogger::info("Update time:   %9u [ms]\n",this->updateStopWatch.getMiliSeconds());

    this->statsCounter++;
}

void RSolverFluidHeat::computeShapeDerivatives()
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

void RSolverFluidHeat::clearShapeDerivatives()
{
    for (uint i=0;i<this->shapeDerivations.size();i++)
    {
        delete this->shapeDerivations[i];
    }
    this->shapeDerivations.clear();
}

void RSolverFluidHeat::computeElement(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager)
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

void RSolverFluidHeat::computeElementGeneral(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();
    uint nInp = RElement::getNIntegrationPoints(element.getType());

    double ro = this->elementDensity[elementID];
    double c = this->elementCapacity[elementID];
    double k = this->elementConduction[elementID];

    double ca = ro*c;

    Ae.fill(0.0);
    be.fill(0.0);

    FluidHeatMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());
    matrixCotainer.clear();

    // Element level matricies
    RRMatrix &me = matrixCotainer.me;
    RRMatrix &ce = matrixCotainer.ce;
    RRMatrix &ke = matrixCotainer.ke;
    RRMatrix &cte = matrixCotainer.cte;
    RRMatrix &kte = matrixCotainer.kte;
    RRMatrix &yte = matrixCotainer.yte;

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

    for (uint intPoint=0;intPoint<nInp;intPoint++)
    {
        const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),intPoint);
        const RRVector &N = shapeFunc.getN();
        const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(intPoint);
        double detJ = this->shapeDerivations[elementID]->getJacobian(intPoint);

        double integValue = detJ * shapeFunc.getW();

        // velocity divergence
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
        // Reynolds numbers
        double Re(k == 0.0 ? 0.0 : mvh * h / (2.0 * k));

        // SUPG stabilization parameter
        double Tsupg = 0.0;
        if (mvh > 0.0)
        {
            Tsupg = h / (2.0 * mvh);
            if (Re > 0.0 && Re <= 3.0)
            {
                Tsupg *= Re / 3.0;
            }
        }

        for (uint m=0;m<nen;m++)
        {
            for (uint n=0;n<nen;n++)
            {
                // m matrix
                if (unsteady)
                {
                    me[m][n] = ca * N[m] * N[n];
                }
                // c matrix
                ce[m][n] = ca * N[m] * vdiv[n];
                // k matrix - the weak form of -div(k*grad(T)), which enters with
                // the same sign as the advection term c.
                ke[m][n] = k * (B[m][0] * B[n][0] + B[m][1] * B[n][1] + B[m][2] * B[n][2]);
                // k~ matrix
                kte[m][n] = Tsupg * ca * vdiv[m] * vdiv[n];
                // y~ matrix
                yte[m][n] = Tsupg * vdiv[m];
            }
            // f vector
            fv[m] = this->nodeHeat[element.getNodeId(m)] * N[m];
        }

        // Assembly element level matrixes
        for (uint m=0;m<nen;m++)
        {
            for (uint n=0;n<nen;n++)
            {
                if (unsteady)
                {
                    Ae[m][n] += me[m][n] + cte[m][n]
                              + alpha * dt * (  ce[m][n]  + ke[m][n]
                                             + kte[m][n] + yte[m][n] );
                }
                else
                {
                    Ae[m][n] += ce[m][n] + ke[m][n]
                              + kte[m][n] + yte[m][n];
                }
            }
            if (unsteady)
            {
                be[m] += dt * fv[m];
                for (uint n=0;n<nen;n++)
                {
                    be[m] += (me[m][n] + cte[m][n] - (1.0 - alpha) * dt * (ce[m][n] + ke[m][n] + kte[m][n] + yte[m][n]))
                           * this->nodeTemperature[element.getNodeId(n)];
                }
            }
            else
            {
                be[m] = fv[m];
            }
        }
        for (uint m=0;m<nen;m++)
        {
            for (uint n=0;n<nen;n++)
            {
                Ae[m][n] *= integValue;
            }
            be[m] *= integValue;
        }
    }
}

void RSolverFluidHeat::computeElementConstantDerivative(unsigned int elementID, RRMatrix &Ae, RRVector &be, RMatrixManager<FluidHeatMatrixContainer> &matrixManager)
{
    bool unsteady = (this->pModel->getTimeSolver().getEnabled());

    const RElement &element = this->pModel->getElement(elementID);
    uint nen = element.size();

    double ro = this->elementDensity[elementID];
    double c = this->elementCapacity[elementID];
    double k = this->elementConduction[elementID];

    double ca = ro*c;

    Ae.fill(0.0);
    be.fill(0.0);

    FluidHeatMatrixContainer &matrixCotainer = matrixManager.getMatricies(element.getType());
    matrixCotainer.clear();

    // Element level matricies
    RRMatrix &me = matrixCotainer.me;
    RRMatrix &ce = matrixCotainer.ce;
    RRMatrix &ke = matrixCotainer.ke;
    RRMatrix &cte = matrixCotainer.cte;
    RRMatrix &kte = matrixCotainer.kte;
    RRMatrix &yte = matrixCotainer.yte;

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

    const RRVector &iN = RElement::getMassVector(element.getType());
    const RRMatrix &iNiN = RElement::getMassMatrix(element.getType());
    double wt = RElement::getTotalWeightFactor(element.getType());
    const RRMatrix &B = this->shapeDerivations[elementID]->getDerivative(0);

    // velocity divergence
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
    // Reynolds numbers
    double Re(k == 0.0 ? 0.0 : mvh * h / (2.0 * k));

    // SUPG stabilization parameter
    double Tsupg = 0.0;
    if (mvh > 0.0)
    {
        Tsupg = h / (2.0 * mvh);
        if (Re > 0.0 && Re <= 3.0)
        {
            Tsupg *= Re / 3.0;
        }
    }

    for (uint m=0;m<nen;m++)
    {
        for (uint n=0;n<nen;n++)
        {
            // m matrix
            if (unsteady)
            {
                me[m][n] = ca * iNiN[m][n];
            }
            // c matrix
            ce[m][n] = ca * iN[m] * vdiv[n];
            // k matrix - the weak form of -div(k*grad(T)), which enters with
            // the same sign as the advection term c.
            ke[m][n] = k * wt * (B[m][0] * B[n][0] + B[m][1] * B[n][1] + B[m][2] * B[n][2]);
            // k~ matrix
            kte[m][n] = Tsupg * ca * vdiv[m] * vdiv[n] * wt;
            // y~ matrix
            yte[m][n] = Tsupg * vdiv[m] * wt;
        }
        // f vector
        fv[m] = this->nodeHeat[element.getNodeId(m)] * iN[m];
    }

    // Assembly element level matrixes
    for (uint m=0;m<nen;m++)
    {
        for (uint n=0;n<nen;n++)
        {
            if (unsteady)
            {
                Ae[m][n] = me[m][n] + cte[m][n]
                         + alpha * dt * (  ce[m][n]  + ke[m][n]
                                         + kte[m][n] + yte[m][n] );
            }
            else
            {
                Ae[m][n] = ce[m][n] + ke[m][n]
                         + kte[m][n] + yte[m][n];
            }
        }
        if (unsteady)
        {
            be[m] = dt * fv[m];
            for (uint n=0;n<nen;n++)
            {
                be[m] += (me[m][n] + cte[m][n] - (1.0 - alpha) * dt * (ce[m][n] + ke[m][n] + kte[m][n] + yte[m][n]))
                       * this->nodeTemperature[element.getNodeId(n)];
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

void RSolverFluidHeat::assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be)
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
                fe[n] -= Ae[n][m] * this->nodeTemperature[rElement.getNodeId(m)];
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

void RSolverFluidHeat::assemblyMatrix(unsigned int elementID, const RRMatrix &Ae, const RRVector &be, RSparseMatrix &Ap, RRVector &bp)
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
                fe[n] -= Ae[n][m] * this->nodeTemperature[nodeID];
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

#include <atomic>
#include <cmath>

#include <omp.h>

#include "rsolverheat.h"
#include "rsolverfluidheat.h"
#include "rconvection.h"
#include "rmatrixsolver.h"

const QString RSolverHeat::solidNodeTemperatureKey("solid-node-temperature");

RSolverHeat::RSolverHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , wallCoupled(false)
    , cvgT(0.0)
{
    this->problemType = R_PROBLEM_HEAT;
}

RSolverHeat::~RSolverHeat()
{
}

bool RSolverHeat::hasConverged() const
{
    // Without walls driven by the fluid heat solver there is nothing to iterate
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

void RSolverHeat::storeSharedData()
{
    this->RSolverGeneric::storeSharedData();

    // Published under a key of its own so the fluid heat solver can hold its
    // walls at the temperature of the solid. Shared in SI units - the fluid heat
    // solver works in scales of its own.
    RRVector temperature(this->nodeTemperature);
    temperature *= 1.0 / this->scales.findScaleFactor(R_VARIABLE_TEMPERATURE);
    this->pSharedData->addData(RSolverHeat::solidNodeTemperatureKey,temperature);
}

void RSolverHeat::recoverSharedData()
{
    this->RSolverGeneric::recoverSharedData();

    // Forced convection walls are driven by the fluid heat solver result.
    this->fluidWallHtc.clear();
    this->fluidWallHtt.clear();
    if (this->pSharedData->hasData(RSolverFluidHeat::wallHeatTransferCoefficientKey,this->pModel->getNElements()) &&
        this->pSharedData->hasData(RSolverFluidHeat::wallFluidTemperatureKey,this->pModel->getNElements()))
    {
        this->fluidWallHtc = this->pSharedData->findData(RSolverFluidHeat::wallHeatTransferCoefficientKey);
        this->fluidWallHtt = this->pSharedData->findData(RSolverFluidHeat::wallFluidTemperatureKey);

        double htcScale = this->scales.findScaleFactor(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
        double temperatureScale = this->scales.findScaleFactor(R_VARIABLE_TEMPERATURE);
        for (uint i=0;i<this->fluidWallHtc.size();i++)
        {
            if (this->fluidWallHtc[i] >= 0.0)
            {
                this->fluidWallHtc[i] *= htcScale;
                this->fluidWallHtt[i] *= temperatureScale;
            }
        }
    }
}

void RSolverHeat::findComputableElements(RProblemType problemType)
{
    this->RSolverGeneric::findComputableElements(problemType);

    // Mark the nodes of solid and of fluid volume elements. A fluid volume is
    // never computed here, whatever its material carries or whichever condition
    // is applied to it - the fluid heat solver owns it.
    RBVector solidNodes(this->pModel->getNNodes(),false);
    RBVector fluidNodes(this->pModel->getNNodes(),false);

    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        const RVolume &rVolume = this->pModel->getVolume(i);
        bool fluid = rVolume.getMaterial().isFluid();

        for (uint j=0;j<rVolume.size();j++)
        {
            uint elementID = rVolume.get(j);
            const RElement &rElement = this->pModel->getElement(elementID);

            if (fluid)
            {
                this->computableElements[elementID] = false;
            }
            else if (!this->computableElements[elementID])
            {
                continue;
            }

            RBVector &nodes = fluid ? fluidNodes : solidNodes;
            for (uint k=0;k<rElement.size();k++)
            {
                nodes[rElement.getNodeId(k)] = true;
            }
        }
    }

    // A point, line or surface entity made computable by a condition, rather
    // than by a solid material of its own, is dropped where it lies inside the
    // fluid - an inlet, an outlet. Nothing connects it to the solid, so its
    // nodes would be left without a stiffness. A wall between the solid and the
    // fluid touches the solid with all its nodes and stays.
    QList<RMaterialProperty::Type> propList;
    for (uint i=0;i<RMaterialProperty::nTypes;i++)
    {
        if (problemType & RMaterialProperty::getProblemTypeMask(RMaterialProperty::Type(i)))
        {
            propList.push_back(RMaterialProperty::Type(i));
        }
    }

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        const RMaterial &rMaterial = pElementGroup->getMaterial();
        if (!rMaterial.isFluid() && rMaterial.hasProperties(propList))
        {
            continue;
        }

        uint nDropped = 0;
        for (uint j=0;j<pElementGroup->size();j++)
        {
            uint elementID = pElementGroup->get(j);
            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &rElement = this->pModel->getElement(elementID);
            if (R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
            {
                continue;
            }
            bool inFluid = true;
            bool onSolid = true;
            for (uint k=0;k<rElement.size();k++)
            {
                inFluid = inFluid && fluidNodes[rElement.getNodeId(k)];
                onSolid = onSolid && solidNodes[rElement.getNodeId(k)];
            }
            if (inFluid && !onSolid)
            {
                this->computableElements[elementID] = false;
                nDropped++;
            }
        }

        if (nDropped == 0)
        {
            continue;
        }

        // Only a condition the fluid heat solver does not read is actually lost.
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            RProblemTypeMask mask = RBoundaryCondition::getProblemTypeMask(pElementGroup->getBoundaryCondition(j).getType());
            if ((mask & R_PROBLEM_HEAT) && !(mask & R_PROBLEM_FLUID_HEAT))
            {
                RLogger::warning("Entity \'%s\' lies inside a fluid domain - its \'%s\' boundary condition is ignored by the heat solver, which solves solids only.\n",
                                 pElementGroup->getName().toUtf8().constData(),
                                 RBoundaryCondition::getName(pElementGroup->getBoundaryCondition(j).getType()).toUtf8().constData());
            }
        }
    }
}

bool RSolverHeat::findFluidWall(uint elementId, double &htc, double &htt) const
{
    if (elementId >= this->fluidWallHtc.size() || elementId >= this->fluidWallHtt.size())
    {
        return false;
    }
    if (this->fluidWallHtc[elementId] < 0.0)
    {
        return false;
    }

    htc = this->fluidWallHtc[elementId];
    htt = this->fluidWallHtt[elementId];

    return true;
}

double RSolverHeat::findTemperatureScale() const
{
    return 1.0;
}

void RSolverHeat::initialize()
{
}

void RSolverHeat::updateScales()
{
    this->scales.setMetre(this->findMeshScale());
    this->scales.setKelvin(this->findTemperatureScale());
}

void RSolverHeat::recover()
{
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
}

void RSolverHeat::prepare()
{
    const bool timeSolverEnabled = this->pModel->getTimeSolver().getEnabled();

    RBVector temperatureSetValues;
    RBVector heatSetValues;
    RBVector heatRateAreaSetValues;
    RBVector heatRateVolumeSetValues;

    // Set again while the surfaces are prepared below.
    this->wallCoupled = false;

    this->generateNodeBook(R_PROBLEM_HEAT);

//    this->pModel->convertNodeToElementVector(this->nodeTemperature,this->elementTemperature);

    this->generateVariableVector(R_VARIABLE_TEMPERATURE,this->elementTemperature,temperatureSetValues,true,this->firstRun,this->firstRun);
    // The Heat boundary condition prescribes the total heat input for the whole
    // entity - it is spread over the entity measure to give the source density
    // the assembly below integrates.
    this->generateHeatVector(this->elementHeat,heatSetValues);
    // Heat rate per unit area / volume are boundary conditions only - they are
    // applied to surface and volume elements respectively, where the source
    // term is already integrated over the element measure.
    this->generateVariableVector(R_VARIABLE_HEAT_RATE_AREA,this->elementHeatRateArea,heatRateAreaSetValues,true,false,false);
    this->generateVariableVector(R_VARIABLE_HEAT_RATE_VOLUME,this->elementHeatRateVolume,heatRateVolumeSetValues,true,false,false);
    this->generateMaterialVecor(RMaterialProperty::ThermalConductivity,this->elementConduction);
    this->generateMaterialVecor(RMaterialProperty::HeatCapacity,this->elementCapacity);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);

    this->pModel->convertElementToNodeVector(this->elementTemperature,temperatureSetValues,this->nodeTemperature,true);

    uint nEnabled = this->nodeBook.getNEnabled();

    this->b.resize(nEnabled);
    this->x.resize(nEnabled);

    this->A.clear();
    this->A.setNRows(nEnabled);
    this->b.fill(0.0);
    this->x.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b once at the end.
    int np = omp_get_max_threads();
    std::vector<RSparseMatrix> Ap(np);
    std::vector<RRVector> bp(np);
    for (int t=0;t<np;t++)
    {
        Ap[t].setNRows(nEnabled);
        bp[t].resize(nEnabled);
        bp[t].fill(0.0);
    }

    // Prepare point elements.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        RPoint &point = this->pModel->getPoint(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(point.size());j++)
        {
            try
            {
                uint elementID = point.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_POINT(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * detJ
                                         * shapeFunc.getW()
                                         * point.getVolume();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            try
            {
                uint elementID = line.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());

                RRMatrix B(element.size(),1);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += dN[m][0]*J[0][0];
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            double kcnd = (B[m][0]*B[n][0]) * line.getCrossArea() * this->elementConduction[elementID];

                            Ke[m][n] += kcnd * detJ * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * detJ
                                         * shapeFunc.getW()
                                         * line.getCrossArea();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

        double surfaceHtc = 0.0;
        double surfaceHtt = 0.0;

        this->getSimpleConvection(surface,surfaceHtc,surfaceHtt);
        this->reportForcedConvection(surface);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            try
            {
                uint elementID = surface.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),2);

                // htc/htt must be per-iteration - the correlated conditions
                // overwrite them per element and the loop runs in parallel.
                double htc = surfaceHtc;
                double htt = surfaceHtt;
                this->getForcedConvection(surface,elementID,htc,htt);
                this->getNaturalConvection(surface,elementID,htc,htt);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]);
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            double kcnd = (B[m][0]*B[n][0]+B[m][1]*B[n][1]) * surface.getThickness() * this->elementConduction[elementID];
                            // Convection
                            double kcnv = N[m] * N[n] * htc;

                            Ke[m][n] += (kcnd + kcnv) * detJ * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * detJ
                                         * shapeFunc.getW()
                                         * surface.getThickness();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementHeatRateArea[elementID] + this->elementRadiativeHeat[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }

                // Convection force
                double elementArea = 0.0;
                if (element.findArea(this->pModel->getNodes(),elementArea))
                {
                    for (unsigned m=0;m<element.size();m++)
                    {
                        fe[m] += htc * htt * elementArea / element.size();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare volume elements.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        RVolume &volume = this->pModel->getVolume(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(volume.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = volume.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                // Conduction
                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            Ke[m][n] += (B[m][0]*B[n][0] + B[m][1]*B[n][1] + B[m][2]*B[n][2])
                                     * this->elementConduction[elementID]
                                     * detJ
                                     * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * detJ
                                         * shapeFunc.getW();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementHeatRateVolume[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Merge per-thread assembly buffers.
    #pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->A.getNRows());i++)
    {
        for (int t=0;t<np;t++)
        {
            this->A.getVector(uint(i)).addVector(Ap[t].getVector(uint(i)));
            this->b[uint(i)] += bp[t][uint(i)];
        }
    }
}

void RSolverHeat::solve()
{
    try
    {
        RLogger::indent();
        RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,1);
        RLogger::unindent();
    }
    catch (const RError &)
    {
        RLogger::unindent();
        throw;
    }

    // Relative size of the change - ||dT|| / ||T|| - so it can be compared
    // against the convergence value of the task group.
    double dtNorm = 0.0;
    double tNorm = 0.0;

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        if (this->nodeBook.getValue(i,position))
        {
            double dt = this->x[position] - this->nodeTemperature[i];
            dtNorm += dt*dt;
            tNorm += this->x[position]*this->x[position];
            this->nodeTemperature[i] = this->x[position];
        }
    }

    this->cvgT = (tNorm > 0.0) ? std::sqrt(dtNorm/tNorm) : 0.0;

    this->pModel->convertNodeToElementVector(this->nodeTemperature,this->elementTemperature);
}

void RSolverHeat::process()
{
//    RLogger::info("Processing results\n");
    double Qx;
    double Qy;
    double Qz;

    // Initialize heat flux vector vector
    this->elementHeatFlux.resize(this->pModel->getNElements(),RR3Vector(0.0,0.0,0.0));

    // Initialize heat transfer coefficient vector. It stays zero on every
    // element which carries no convection boundary condition.
    this->elementHeatTransferCoefficient.resize(this->pModel->getNElements());
    this->elementHeatTransferCoefficient.fill(0.0);

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        for (uint j=0;j<line.size();j++)
        {
            uint elementID = line.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRVector B(element.size());

            Qx = Qy = Qz = 0.0;

            // Conduction
            if (line.getCrossArea() > 0.0)
            {
                double Qi = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    // Heat flux is a density - the shape function derivatives
                    // are averaged over the integration points and must not be
                    // weighted by the Jacobian determinant.
                    this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    if (line.getCrossArea() != 0.0)
                    {
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m] += dN[m][0] * J[0][0] / double(nInp);
                        }
                    }
                }

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    Qi -= B[k] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                Qx += R[0][0]*Qi;
                Qy += R[1][0]*Qi;
                Qz += R[2][0]*Qi;
            }

            this->elementHeatFlux[elementID][0] = Qx;
            this->elementHeatFlux[elementID][1] = Qy;
            this->elementHeatFlux[elementID][2] = Qz;
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

        double surfaceHtc = 0.0;
        double surfaceHtt = 0.0;

        this->getSimpleConvection(surface,surfaceHtc,surfaceHtt);

        for (uint j=0;j<surface.size();j++)
        {
            uint elementID = surface.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRMatrix B(element.size(),2);

            double htc = surfaceHtc;
            double htt = surfaceHtt;
            this->getForcedConvection(surface,elementID,htc,htt);
            this->getNaturalConvection(surface,elementID,htc,htt);

            this->elementHeatTransferCoefficient[elementID] = htc;

            Qx = Qy = Qz = 0.0;

            // Conduction
            if (surface.getThickness() > 0.0)
            {
                double Qi = 0.0;
                double Qj = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    // Heat flux is a density - see the line element loop above.
                    this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    if (surface.getThickness() != 0.0)
                    {
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]) / double(nInp);
                            B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]) / double(nInp);
                        }
                    }
                }

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    Qi -= B[k][0] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                    Qj -= B[k][1] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                Qx += R[0][0]*Qi + R[0][1]*Qj;
                Qy += R[1][0]*Qi + R[1][1]*Qj;
                Qz += R[2][0]*Qi + R[2][1]*Qj;
            }

            // Convection - Newton's law of cooling gives a flux density, so the
            // element area must not enter here.
            if (htc > 0.0)
            {
                RR3Vector normal;
                if (element.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]))
                {
                    double Qhe = htc * (htt - this->elementTemperature[elementID]);
                    Qx += normal[0] * Qhe;
                    Qy += normal[1] * Qhe;
                    Qz += normal[2] * Qhe;
                }
            }


            this->elementHeatFlux[elementID][0] = Qx;
            this->elementHeatFlux[elementID][1] = Qy;
            this->elementHeatFlux[elementID][2] = Qz;
        }
    }

    // Process volume elements.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        RVolume &volume = this->pModel->getVolume(i);

        for (uint j=0;j<volume.size();j++)
        {
            uint elementID = volume.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
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
                // Heat flux is a density - see the line element loop above.
                this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                for (uint m=0;m<dN.getNRows();m++)
                {
                    B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]) / double(nInp);
                    B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]) / double(nInp);
                    B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]) / double(nInp);
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
    }
}

void RSolverHeat::store()
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

    // Heat transfer coefficient
    uint heatTransferCoefficientPos = this->pModel->findVariable(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
    if (heatTransferCoefficientPos == RConstants::eod)
    {
        heatTransferCoefficientPos = this->pModel->addVariable(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
        this->pModel->getVariable(heatTransferCoefficientPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementHeatTransferCoefficient),
                    RStatistics::findMaximumValue(this->elementHeatTransferCoefficient));
    }
    RVariable &heatTransferCoefficient = this->pModel->getVariable(heatTransferCoefficientPos);

    heatTransferCoefficient.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    heatTransferCoefficient.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        heatTransferCoefficient.setValue(0,i,this->elementHeatTransferCoefficient[i]);
    }

    RLogger::unindent();
}

void RSolverHeat::statistics()
{
    this->printStats(R_VARIABLE_TEMPERATURE);
    this->printStats(R_VARIABLE_HEAT_FLUX);
    this->printStats(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
    this->processMonitoringPoints();

    if (this->wallCoupled)
    {
        RLogger::info("Walls driven by the fluid heat solver - convergence-T: %-13g\n",this->cvgT);
        if (this->taskCvgValue > 0.0)
        {
            RLogger::info("Convergence target: %-13g%s\n",
                          this->taskCvgValue,
                          this->hasConverged() ? " (reached)" : "");
        }
    }
}

void RSolverHeat::assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    const RElement &element = this->pModel->getElement(elementID);

    RRMatrix Ae(element.size(),element.size());
    RRVector be(element.size());

    Ae.fill(0.0);
    be.fill(0.0);

    if (this->pModel->getTimeSolver().getEnabled())
    {
        for (unsigned m=0;m<element.size();m++)
        {
            be[m] = dt * fe[m];
            for (unsigned n=0;n<element.size();n++)
            {
                Ae[m][n] = Me[m][n] + alpha * dt * Ke[m][n];
                be[m] += (Me[m][n] - (1.0 - alpha) * dt * Ke[m][n]) * this->nodeTemperature[element.getNodeId(n)];
            }
        }
    }
    else
    {
        Ae = Ke;
        be = fe;
    }

    // Apply explicit boundary conditions.
    for (uint m=0;m<element.size();m++)
    {
        uint position;
        uint nodeID = element.getNodeId(m);
        if (!this->nodeBook.getValue(nodeID,position))
        {
            for (uint n=0;n<element.size();n++)
            {
                be[n] -= Ae[n][m] * this->nodeTemperature[nodeID];
            }
        }
    }


    // Assembly final matrix system
    for (uint m=0;m<element.size();m++)
    {
        uint mp;

        if (this->nodeBook.getValue(element.getNodeId(m),mp))
        {
            bp[mp] += be[m];
            for (uint n=0;n<element.size();n++)
            {
                uint np = 0;

                if (this->nodeBook.getValue(element.getNodeId(n),np))
                {
                    Ap.addValue(mp,np,Ae[m][n]);
                }
            }
        }
    }
}

bool RSolverHeat::getSimpleConvection(const RElementGroup &elementGroup, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE))
    {
        return false;
    }

    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE);
    uint cPos = 0;

    cPos = bc.findComponentPosition(R_VARIABLE_CONVECTION_COEFFICIENT);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_CONVECTION_COEFFICIENT).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE).toUtf8().constData());
    }
    htc = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    cPos = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_FLUID_TEMPERATURE).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE).toUtf8().constData());
    }
    htt = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    return true;
}

void RSolverHeat::checkConvectionInput(double value,
                                       RVariableType variableType,
                                       RBoundaryConditionType boundaryConditionType,
                                       const RElementGroup &elementGroup) const
{
    if (std::fabs(value) > RConstants::eps)
    {
        return;
    }

    throw RError(RError::Type::Application,R_ERROR_REF,
                 "Value of \'%s\' configured in \'%s\' boundary condition on entity \'%s\' is zero - the convection correlation can not be evaluated.",
                 RVariable::getName(variableType).toUtf8().constData(),
                 RBoundaryCondition::getName(boundaryConditionType).toUtf8().constData(),
                 elementGroup.getName().toUtf8().constData());
}

void RSolverHeat::reportForcedConvection(const RElementGroup &elementGroup)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED))
    {
        return;
    }

    uint nFromFluid = 0;
    uint nCorrelated = 0;
    uint nComputable = 0;
    bool correlationChecked = false;
    bool correlationApplies = false;

    for (uint i=0;i<elementGroup.size();i++)
    {
        uint elementID = elementGroup.get(i);
        if (!this->computableElements[elementID])
        {
            continue;
        }
        nComputable++;

        double htc = 0.0;
        double htt = 0.0;
        if (this->findFluidWall(elementID,htc,htt))
        {
            nFromFluid++;
        }
        else
        {
            // The configured values are the same on every element, so a single
            // evaluation tells whether the correlation applies at all - and
            // stops on a value it can not work with.
            if (!correlationChecked)
            {
                correlationApplies = this->getForcedConvection(elementGroup,elementID,htc,htt);
                correlationChecked = true;
            }
            if (correlationApplies)
            {
                nCorrelated++;
            }
        }
    }

    if (nFromFluid > 0)
    {
        this->wallCoupled = true;
    }

    if (nComputable == 0)
    {
        return;
    }
    if (nFromFluid == nComputable)
    {
        RLogger::info("Forced convection on entity \'%s\' takes the heat transfer computed by the fluid heat solver.\n",
                      elementGroup.getName().toUtf8().constData());
        return;
    }
    if (nFromFluid + nCorrelated == 0)
    {
        RLogger::warning("Forced convection on entity \'%s\' is ignored - no fluid heat result covers this surface and the condition carries no fluid temperature to fall back on.\n",
                         elementGroup.getName().toUtf8().constData());
        return;
    }
    if (nFromFluid == 0)
    {
        RLogger::info("Forced convection on entity \'%s\' is correlated from the values configured on the condition - no fluid heat result covers this surface.\n",
                      elementGroup.getName().toUtf8().constData());
        return;
    }
    RLogger::info("Forced convection on entity \'%s\' takes the heat transfer computed by the fluid heat solver on %u of %u elements, the rest is correlated from the values configured on the condition.\n",
                  elementGroup.getName().toUtf8().constData(),
                  nFromFluid,
                  nComputable);
}

bool RSolverHeat::getForcedConvection(const RElementGroup &elementGroup, uint elementId, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED))
    {
        return false;
    }

    // A wall bordering a meshed fluid takes the heat transfer the fluid heat
    // solver computed behind it. No correlation is involved - the resolved flow
    // already carries the convection.
    if (this->findFluidWall(elementId,htc,htt))
    {
        return true;
    }

    // The correlation is the fall-back for a surface no such result covers - one
    // bordering no meshed fluid, or a run where the fluid heat solver has not
    // produced a result yet.
    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED);
    uint cPos = 0;

    uint fluidTemperaturePosition = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
    if (fluidTemperaturePosition == RConstants::eod)
    {
        // Conditions stored before the component existed carry no value to
        // fall back on, so there is nothing to exchange heat with.
        return false;
    }
    htt = bc.getComponent(fluidTemperaturePosition).get(this->pModel->getTimeSolver().getCurrentTime());

    // Density
    cPos = bc.findComponentPosition(R_VARIABLE_DENSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DENSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double ro = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(ro,R_VARIABLE_DENSITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Dynamic viscosity
    cPos = bc.findComponentPosition(R_VARIABLE_DYNAMIC_VISCOSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DYNAMIC_VISCOSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double mu = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(mu,R_VARIABLE_DYNAMIC_VISCOSITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Heat capacity
    cPos = bc.findComponentPosition(R_VARIABLE_HEAT_CAPACITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HEAT_CAPACITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double c = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(c,R_VARIABLE_HEAT_CAPACITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Hydraulic diameter
    cPos = bc.findComponentPosition(R_VARIABLE_HYDRAULIC_DIAMETER);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HYDRAULIC_DIAMETER).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double d = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(d,R_VARIABLE_HYDRAULIC_DIAMETER,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Thermal conductivity
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_CONDUCTIVITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_CONDUCTIVITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double k = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(k,R_VARIABLE_THERMAL_CONDUCTIVITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Mean velocity
    cPos = bc.findComponentPosition(R_VARIABLE_VELOCITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_VELOCITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double v = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(v,R_VARIABLE_VELOCITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    RConvection convection;

    convection.setType(R_CONVECTION_FORCED_EXTERNAL);
    convection.setMaterial("Custom material",mu,ro,k,c,0.0);
    convection.setDiameter(d);
    convection.setVelocity(v);
    convection.setFluidTemp(htt);

    htc = convection.calculateHtc();

    return true;
}

bool RSolverHeat::getNaturalConvection(const RElementGroup &elementGroup, uint elementId, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_NATURAL))
    {
        return false;
    }

    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_NATURAL);
    uint cPos = 0;

    // Density
    cPos = bc.findComponentPosition(R_VARIABLE_DENSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DENSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double ro = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    // Dynamic viscosity
    cPos = bc.findComponentPosition(R_VARIABLE_DYNAMIC_VISCOSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DYNAMIC_VISCOSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double mu = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(mu,R_VARIABLE_DYNAMIC_VISCOSITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Fluid temperature
    cPos = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_FLUID_TEMPERATURE).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    htt = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    // Heat capacity
    cPos = bc.findComponentPosition(R_VARIABLE_HEAT_CAPACITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HEAT_CAPACITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double c = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(c,R_VARIABLE_HEAT_CAPACITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Hydraulic diameter
    cPos = bc.findComponentPosition(R_VARIABLE_HYDRAULIC_DIAMETER);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HYDRAULIC_DIAMETER).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double d = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(d,R_VARIABLE_HYDRAULIC_DIAMETER,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Thermal conductivity
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_CONDUCTIVITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_CONDUCTIVITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double k = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(k,R_VARIABLE_THERMAL_CONDUCTIVITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Thermal expansion coefficient
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_EXPANSION_COEFFICIENT);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_EXPANSION_COEFFICIENT).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double b = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    RConvection convection;

    convection.setType(R_CONVECTION_NATURAL_EXTERNAL_HORIZONTAL_PLATES);
    convection.setMaterial("Custom material",mu,ro,k,c,b);
    convection.setDiameter(d);
    convection.setSurfTemp(this->elementTemperature[elementId]);
    convection.setFluidTemp(htt);

    htc = convection.calculateHtc();

    return true;
}

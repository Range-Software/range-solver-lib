#include <algorithm>
#include <atomic>
#include <cmath>

#include <omp.h>

#include "rsolverstress.h"
#include "rmatrixsolver.h"
#include "reigenvaluesolver.h"

namespace
{

//! RElement::findJacobian() returns for a surface element a 3n x 2n matrix with
//! the element frame in every node block. Rotating an element matrix between
//! the element and the global frame needs the block diagonal part only - the
//! off-diagonal blocks would add up the couplings of all node pairs.
RRMatrix findSurfaceRotation(const RRMatrix &Rt)
{
    uint nn = Rt.getNRows() / 3;
    RRMatrix R(3*nn,2*nn,0.0);
    for (uint i=0;i<nn;i++)
    {
        for (uint r=0;r<3;r++)
        {
            R[3*i+r][2*i+0] = Rt[3*i+r][2*i+0];
            R[3*i+r][2*i+1] = Rt[3*i+r][2*i+1];
        }
    }
    return R;
}

}

RSolverStress::RSolverStress(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData, bool modalAnalysis)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , nodeLoadOldValid(false)
    , oldStateTimeStep(RConstants::eod)
    , dynamic(false)
    , computeInitialAcceleration(false)
    , newmarkGamma(0.5)
    , newmarkBeta(0.25)
    , hhtAlpha(0.0)
    , timeStepSize(0.0)
{
    this->problemType = modalAnalysis ? R_PROBLEM_STRESS_MODAL : R_PROBLEM_STRESS;
}

RSolverStress::~RSolverStress()
{

}

bool RSolverStress::hasConverged() const
{
    return true;
}

uint RSolverStress::getNComputedModes() const
{
    return uint(this->d.size());
}

void RSolverStress::initialize()
{
}

void RSolverStress::updateScales()
{
    this->scales.setMetre(this->findMeshScale());
}

void RSolverStress::recover()
{
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeDisplacement.x,0.0);
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeDisplacement.y,0.0);
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeDisplacement.z,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeForce.x,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeForce.y,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeForce.z,0.0);
    this->recoverVariable(R_VARIABLE_PRESSURE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodePressure,0.0);
    // Velocity and acceleration carry the Newmark state over to the next time
    // step and over a restart.
    this->recoverVariable(R_VARIABLE_STRUCTURAL_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
    this->recoverVariable(R_VARIABLE_STRUCTURAL_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
    this->recoverVariable(R_VARIABLE_STRUCTURAL_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
    this->recoverVariable(R_VARIABLE_STRUCTURAL_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeAcceleration.x,0.0);
    this->recoverVariable(R_VARIABLE_STRUCTURAL_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeAcceleration.y,0.0);
    this->recoverVariable(R_VARIABLE_STRUCTURAL_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeAcceleration.z,0.0);

//    this->syncShared("node-displacement-x",this->nodeDisplacement.x);
//    this->syncShared("node-displacement-y",this->nodeDisplacement.y);
//    this->syncShared("node-displacement-z",this->nodeDisplacement.z);

}

void RSolverStress::prepare()
{
    this->dynamic = this->findDynamic();
    this->computeInitialAcceleration = false;

    const bool needsMass = this->dynamic || this->problemType == R_PROBLEM_STRESS_MODAL;

    //! Element force vector.
    RSolverCartesianVector<RRVector> elementForce;
    RSolverCartesianVector<RBVector> forceSetValues;


    //! Element gravity vector.
    RSolverCartesianVector<RRVector> elementGravity;
    RSolverCartesianVector<RBVector> gGravitySetValues;

    //! Element pressure.
    RRVector elementPressure;
    RBVector pressureSetValues;

    //! Element traction per unit area.
    RSolverCartesianVector<RRVector> elementForceUnitArea;
    RSolverCartesianVector<RBVector> forceUnitAreaSetValues;

    //! Element weight.
    RRVector elementWeight;
    RBVector weightSetValues;

    RBVector temperatureSetValues;

    // The node frames and the directions they hold come from the constraints
    // themselves, so they have to be known before the node book is built.
    this->generateLocalConstraints();
    this->generateNodeBook();
    this->generateVariableVector(R_VARIABLE_FORCE_X,elementForce.x,forceSetValues.x,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_FORCE_Y,elementForce.y,forceSetValues.y,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_FORCE_Z,elementForce.z,forceSetValues.z,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_X,elementGravity.x,gGravitySetValues.x,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Y,elementGravity.y,gGravitySetValues.y,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Z,elementGravity.z,gGravitySetValues.z,true,true,true);
    this->generateVariableVector(R_VARIABLE_PRESSURE,elementPressure,pressureSetValues,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_X,elementForceUnitArea.x,forceUnitAreaSetValues.x,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_Y,elementForceUnitArea.y,forceUnitAreaSetValues.y,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_Z,elementForceUnitArea.z,forceUnitAreaSetValues.z,true,true,true);
    this->generateVariableVector(R_VARIABLE_WEIGHT,elementWeight,weightSetValues,true,true,true);
    this->generateMaterialVecor(RMaterialProperty::ModulusOfElasiticity,this->elementElasticityModulus);
    this->generateMaterialVecor(RMaterialProperty::PoissonRatio,this->elementPoissonRatio);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);
    this->generateMaterialVecor(RMaterialProperty::ThermalExpansionCoefficient,this->elementThermalExpansion);
    this->generateVariableVector(R_VARIABLE_TEMPERATURE,this->elementEnvironmentTemperature,temperatureSetValues,false,false,true);
    this->generateMaterialVecor(RMaterialProperty::RayleighMassDamping,this->elementRayleighMass);
    this->generateMaterialVecor(RMaterialProperty::RayleighStiffnessDamping,this->elementRayleighStiffness);

    if (this->dynamic)
    {
        this->findTimeIntegrationParameters();
        if (this->firstRun)
        {
            this->applyInitialConditions();
        }
        this->storeOldState();
    }

    uint nn = this->pModel->getNNodes();

    this->b.resize(this->nodeBook.getNEnabled());
    this->x.resize(this->nodeBook.getNEnabled());

    this->M.clear();
    this->M.setNRows(this->b.size());
    this->D.clear();
    this->D.setNRows(this->computeInitialAcceleration ? this->b.size() : 0);
    this->A.clear();
    this->A.setNRows(this->b.size());
    this->b.fill(0.0);
    this->x.fill(0.0);
    this->r.resize(this->computeInitialAcceleration ? this->b.size() : 0);
    this->r.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b (and M for modal) at the end.
    int np = omp_get_max_threads();
    std::vector<AssemblyBuffer> buffers(np);
    for (int t=0;t<np;t++)
    {
        buffers[t].A.setNRows(this->b.size());
        buffers[t].M.setNRows(this->b.size());
        buffers[t].D.setNRows(this->D.getNRows());
        buffers[t].b.resize(this->b.size(),0.0);
        buffers[t].b.fill(0.0);
        buffers[t].r.resize(this->r.size(),0.0);
        buffers[t].r.fill(0.0);
        buffers[t].f.resize(this->dynamic ? 3*nn : 0,0.0);
        buffers[t].f.fill(0.0);
    }

    this->pModel->convertElementToNodeVector(elementForce.x,forceSetValues.x,this->nodeForce.x,true);
    this->pModel->convertElementToNodeVector(elementForce.y,forceSetValues.y,this->nodeForce.y,true);
    this->pModel->convertElementToNodeVector(elementForce.z,forceSetValues.z,this->nodeForce.z,true);
    this->pModel->convertElementToNodeVector(elementPressure,pressureSetValues,this->nodePressure,true);

    // Convert node pressure to element pressure
    for (uint i=0;i<elementPressure.size();i++)
    {
        if (!pressureSetValues[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            uint nne = rElement.size();
            if (nne > 0)
            {
                elementPressure[i] = 0.0;
                for (uint j=0;j<nne;j++)
                {
                    elementPressure[i] += this->nodePressure[rElement.getNodeId(j)];
                }
                elementPressure[i] /= double(nne);
            }
        }
    }

    std::vector<RNode> nodesBkp;
    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        this->nodeInitialDisplacement.x = this->nodeDisplacement.x;
        this->nodeInitialDisplacement.y = this->nodeDisplacement.y;
        this->nodeInitialDisplacement.z = this->nodeDisplacement.z;
        RLogger::info("Moving prestressed nodes\n");
        nodesBkp = this->pModel->getNodes();
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            this->pModel->getNode(i).move(RR3Vector(this->nodeDisplacement.x[i],this->nodeDisplacement.y[i],this->nodeDisplacement.z[i]));
        }
    }

    // Prepare point elements.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        RPoint &point = this->pModel->getPoint(i);
        double pointVolume = point.getVolume();
        // Force and Weight are totals over the entity, exactly as they are for
        // a line or a surface, so they are spread over its point elements.
        double pointCount = double(std::max(point.size(),uint(1)));

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(point.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = point.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_POINT(this->pModel->getElement(elementID).getType()));
                RRMatrix Me(3,3);
                RRMatrix Ke(3,3);
                RRVector fe(3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                // Force
                fe[0] += elementForce.x[elementID] / pointCount;
                fe[1] += elementForce.y[elementID] / pointCount;
                fe[2] += elementForce.z[elementID] / pointCount;
                // Weight
                fe[0] += elementWeight[elementID] * elementGravity.x[elementID] / pointCount;
                fe[1] += elementWeight[elementID] * elementGravity.y[elementID] / pointCount;
                fe[2] += elementWeight[elementID] * elementGravity.z[elementID] / pointCount;
                // Own weight
                if (pointVolume > 0.0)
                {
                    fe[0] += elementGravity.x[elementID] * this->elementDensity[elementID] * pointVolume;
                    fe[1] += elementGravity.y[elementID] * this->elementDensity[elementID] * pointVolume;
                    fe[2] += elementGravity.z[elementID] * this->elementDensity[elementID] * pointVolume;
                }

                // Mass
                if (needsMass)
                {
                    Me.setIdentity(3);
                    Me *= this->elementDensity[elementID] * pointVolume;
                }

                this->assemblyMatrix(elementID,Me,Ke,fe,buffers[uint(omp_get_thread_num())]);
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
        double lineCrossArea = line.getCrossArea();

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = line.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3,0.0);
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);

                RRMatrix Be(3*element.size(),1);
                RRMatrix BeT(1,3*element.size());

                double lineLength = 0.0;
                element.findLength(this->pModel->getNodes(),lineLength);

                double E = this->elementElasticityModulus[elementID];
                double De = E * lineCrossArea;

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    if (lineCrossArea > 0.0)
                    {
                        // Strain-displacement vector of the truss - it maps the
                        // global nodal displacements onto the axial strain.
                        Be.fill(0.0);
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            Be[3*m+0][0] = Rt[3*m+0][0]*dN[m][0]*J[0][0];
                            Be[3*m+1][0] = Rt[3*m+1][0]*dN[m][0]*J[0][0];
                            Be[3*m+2][0] = Rt[3*m+2][0]*dN[m][0]*J[0][0];
                        }
                        BeT.transpose(Be);

                        // Ke += E*A * B * B^T * detJ * W, accumulated over the
                        // integration points.
                        RRMatrix BeScaled(Be);
                        BeScaled *= De * detJ * shapeFunc.getW();
                        RRMatrix::mlt(BeScaled,BeT,Ke,true);

                        // Thermal expansion force: f += E*A * alpha * dT * B * detJ * W
                        double thermalFactor = De
                                             * this->elementThermalExpansion[elementID]
                                             * dT
                                             * detJ
                                             * shapeFunc.getW();
                        for (uint m=0;m<3*element.size();m++)
                        {
                            fe[m] += thermalFactor * Be[m][0];
                        }
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        if (lineCrossArea > 0.0)
                        {
                            // Mass
                            if (needsMass)
                            {
                                for (uint n=0;n<element.size();n++)
                                {
                                    double value = N[m] * N[n]
                                                 * this->elementDensity[elementID]
                                                 * detJ
                                                 * shapeFunc.getW()
                                                 * lineCrossArea;
                                    Me[3*m+0][3*n+0] += std::pow(Rt[0][0],2.0)*value;
                                    Me[3*m+1][3*n+1] += std::pow(Rt[1][0],2.0)*value;
                                    Me[3*m+2][3*n+2] += std::pow(Rt[2][0],2.0)*value;
                                }
                            }
                        }

                        double integValue = N[m] * detJ * shapeFunc.getW();

                        // Force
                        fe[3*m+0] += (elementForce.x[elementID] / lineLength) * integValue;
                        fe[3*m+1] += (elementForce.y[elementID] / lineLength) * integValue;
                        fe[3*m+2] += (elementForce.z[elementID] / lineLength) * integValue;
                        // Weight
                        fe[3*m+0] += (elementWeight[elementID] * elementGravity.x[elementID] / lineLength) * integValue;
                        fe[3*m+1] += (elementWeight[elementID] * elementGravity.y[elementID] / lineLength) * integValue;
                        fe[3*m+2] += (elementWeight[elementID] * elementGravity.z[elementID] / lineLength) * integValue;
                        // Own weight
                        if (lineCrossArea > 0.0)
                        {
                            fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                            fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                            fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                        }

                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,buffers[uint(omp_get_thread_num())]);
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
        double surfaceArea = surface.findArea(this->pModel->getNodes(),this->pModel->getElements());
        double surfaceThickness = surface.getThickness();

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = surface.get(uint(j));

                if (!this->computableElements[elementID] && !this->includableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3);
                RRMatrix Ke(element.size()*3,element.size()*3);
                RRVector fe(element.size()*3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*2,3);
                RRMatrix BeT(3,element.size()*2);
                RRMatrix BeD(element.size()*2,3);
                RRMatrix Met(element.size()*2,element.size()*2);
                RRMatrix MeRt(element.size()*3,element.size()*2);
                RRMatrix Ket(element.size()*2,element.size()*2);
                RRMatrix KeRt(element.size()*3,element.size()*2);
                RRVector fet(element.size()*2);

                RRMatrix De(3,3,0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1-v;   De[0][1] = v;
                De[1][0] = v;     De[1][1] = 1-v;
                De[2][2] = (1-2*v)/2;
                De *= E/((1+v)*(1-2*v));

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                RR3Vector normal;
                element.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    Rt = findSurfaceRotation(Rt);
                    RtT.transpose(Rt);

                    // Mass of this integration point only - it is rotated to the
                    // global frame and added to Me below.
                    Met.fill(0.0);

                    if (surfaceThickness > 0.0)
                    {
                        B.fill(0.0);
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]);
                            B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]);
                        }

                        for (uint m=0;m<element.size();m++)
                        {
                            Be[2*m][0] = B[m][0];   Be[2*m+1][0] = 0.0;
                            Be[2*m][1] = 0.0;       Be[2*m+1][1] = B[m][1];
                            Be[2*m][2] = B[m][1];   Be[2*m+1][2] = B[m][0];
                        }
                        BeT.transpose(Be);

                        // Ke += Rt * B * D * B^T * Rt^T * t * detJ * W,
                        // accumulated over the integration points. The thickness
                        // enters as it does for the mass and the thermal load.
                        RRMatrix::mlt(Be,De,BeD);
                        RRMatrix::mlt(BeD,BeT,Ket);
                        Ket *= surfaceThickness * detJ * shapeFunc.getW();
                        RRMatrix::mlt(Rt,Ket,KeRt);
                        RRMatrix::mlt(KeRt,RtT,Ke,true);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        if (surfaceThickness > 0.0)
                        {
                            // Mass
                            if (needsMass)
                            {
                                for (uint n=0;n<element.size();n++)
                                {
                                    double value = N[m] * N[n]
                                                 * this->elementDensity[elementID]
                                                 * detJ
                                                 * shapeFunc.getW()
                                                 * surfaceThickness;
                                    Met[2*m+0][2*n+0] += value;
                                    Met[2*m+1][2*n+1] += value;
                                }
                            }
                        }

                        double integValue = N[m] * detJ * shapeFunc.getW();

                        // Pressure vector
                        fe[3*m+0] += elementPressure[elementID] * normal[0] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                        fe[3*m+1] += elementPressure[elementID] * normal[1] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                        fe[3*m+2] += elementPressure[elementID] * normal[2] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                        // Force per unit area
                        fe[3*m+0] += elementForceUnitArea.x[elementID] * integValue;
                        fe[3*m+1] += elementForceUnitArea.y[elementID] * integValue;
                        fe[3*m+2] += elementForceUnitArea.z[elementID] * integValue;
                        // Force
                        fe[3*m+0] += (elementForce.x[elementID] / surfaceArea) * integValue;
                        fe[3*m+1] += (elementForce.y[elementID] / surfaceArea) * integValue;
                        fe[3*m+2] += (elementForce.z[elementID] / surfaceArea) * integValue;
                        // Weight
                        fe[3*m+0] += (elementWeight[elementID] * elementGravity.x[elementID] / surfaceArea) * integValue;
                        fe[3*m+1] += (elementWeight[elementID] * elementGravity.y[elementID] / surfaceArea) * integValue;
                        fe[3*m+2] += (elementWeight[elementID] * elementGravity.z[elementID] / surfaceArea) * integValue;
                        // Own weight
                        if (surfaceThickness > 0.0)
                        {
                            fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * surfaceThickness * integValue;
                            fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * surfaceThickness * integValue;
                            fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * surfaceThickness * integValue;
                        }

                        // Thermal expansion - f = B^T * D * eps_th. In plane
                        // strain the free out of plane expansion is suppressed,
                        // which raises the in-plane thermal strain to
                        // (1+v)*alpha*dT, and there is no thermal shear strain.
                        if (surfaceThickness > 0.0)
                        {
                            double thermalStrain = (1.0 + v) * this->elementThermalExpansion[elementID] * dT;
                            fet.fill(0.0);
                            for (uint n=0;n<2;n++)
                            {
                                fet[2*m+0] += thermalStrain * BeD[2*m+0][n] * surfaceThickness * detJ * shapeFunc.getW();
                                fet[2*m+1] += thermalStrain * BeD[2*m+1][n] * surfaceThickness * detJ * shapeFunc.getW();
                            }

                            fe[3*m+0] += Rt[3*m+0][2*m+0]*fet[2*m+0] + Rt[3*m+0][2*m+1]*fet[2*m+1];
                            fe[3*m+1] += Rt[3*m+1][2*m+0]*fet[2*m+0] + Rt[3*m+1][2*m+1]*fet[2*m+1];
                            fe[3*m+2] += Rt[3*m+2][2*m+0]*fet[2*m+0] + Rt[3*m+2][2*m+1]*fet[2*m+1];
                        }
                    }

                    // Mass
                    if (surfaceThickness > 0.0 && needsMass)
                    {
                        RRMatrix::mlt(Rt,Met,MeRt);
                        RRMatrix::mlt(MeRt,RtT,Me,true);
                    }
                    if (!this->computableElements[elementID])
                    {
                        Me.fill(0.0);
                        Ke.fill(0.0);
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,buffers[uint(omp_get_thread_num())]);
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
                uint elementID = volume.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3);
                RRMatrix Ke(element.size()*3,element.size()*3);
                RRVector fe(element.size()*3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*3,6);
                RRMatrix BeT(6,element.size()*3);
                RRMatrix BeD(element.size()*3,6);
                RRMatrix Ket(element.size()*3,element.size()*3);

                RRMatrix De(6,6);
                De.fill(0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1.0-v; De[0][1] = v;     De[0][2] = v;
                De[1][0] = v;     De[1][1] = 1.0-v; De[1][2] = v;
                De[2][0] = v;     De[2][1] = v;     De[2][2] = 1.0-v;
                De[3][3] = De[4][4] = De[5][5] = (1.0-2.0*v)/2.0;
                De *= E/((1.0+v)*(1.0-2.0*v));

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        Be[3*m+0][0] = B[m][0];   Be[3*m+1][0] = 0.0;       Be[3*m+2][0] = 0.0;
                        Be[3*m+0][1] = 0.0;       Be[3*m+1][1] = B[m][1];   Be[3*m+2][1] = 0.0;
                        Be[3*m+0][2] = 0.0;       Be[3*m+1][2] = 0.0;       Be[3*m+2][2] = B[m][2];
                        Be[3*m+0][3] = 0.0;       Be[3*m+1][3] = B[m][2];   Be[3*m+2][3] = B[m][1];
                        Be[3*m+0][4] = B[m][2];   Be[3*m+1][4] = 0.0;       Be[3*m+2][4] = B[m][0];
                        Be[3*m+0][5] = B[m][1];   Be[3*m+1][5] = B[m][0];   Be[3*m+2][5] = 0.0;
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    for (uint m=0;m<3*element.size();m++)
                    {
                        for (uint n=0;n<3*element.size();n++)
                        {
                            // Stiffness matrix
                            Ke[m][n] += Ket[m][n] * detJ * shapeFunc.getW();
                        }
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<element.size();n++)
                        {
                            // Mass
                            if (needsMass)
                            {
                                double value = N[m] * N[n]
                                             * this->elementDensity[elementID]
                                             * detJ
                                             * shapeFunc.getW();
                                Me[3*m+0][3*n+0] += value;
                                Me[3*m+1][3*n+1] += value;
                                Me[3*m+2][3*n+2] += value;
                            }
                        }

                        // Own weight
                        fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();
                        fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();
                        fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();

                        // Thermal expansion
                        for (uint n=0;n<3;n++)
                        {
                            fe[3*m+0] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+0][n] * detJ * shapeFunc.getW();
                            fe[3*m+1] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+1][n] * detJ * shapeFunc.getW();
                            fe[3*m+2] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+2][n] * detJ * shapeFunc.getW();
                        }
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,buffers[uint(omp_get_thread_num())]);
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
            this->A.getVector(uint(i)).addVector(buffers[t].A.getVector(uint(i)));
            this->M.getVector(uint(i)).addVector(buffers[t].M.getVector(uint(i)));
            this->b[uint(i)] += buffers[t].b[uint(i)];
            if (this->computeInitialAcceleration)
            {
                this->D.getVector(uint(i)).addVector(buffers[t].D.getVector(uint(i)));
                this->r[uint(i)] += buffers[t].r[uint(i)];
            }
        }
    }

    if (this->dynamic)
    {
        this->nodeLoad.x.resize(nn,0.0);
        this->nodeLoad.y.resize(nn,0.0);
        this->nodeLoad.z.resize(nn,0.0);
        this->nodeLoad.x.fill(0.0);
        this->nodeLoad.y.fill(0.0);
        this->nodeLoad.z.fill(0.0);
        for (int t=0;t<np;t++)
        {
            for (uint i=0;i<nn;i++)
            {
                this->nodeLoad.x[i] += buffers[t].f[3*i+0];
                this->nodeLoad.y[i] += buffers[t].f[3*i+1];
                this->nodeLoad.z[i] += buffers[t].f[3*i+2];
            }
        }

        // HHT load term -alpha*f(n). The load of the current time step enters
        // the element vectors as (1+alpha)*f(n+1). Without the load of the
        // previous step (first step, restart) the load is taken as constant
        // over the step.
        if (this->hhtAlpha != 0.0)
        {
            const RSolverCartesianVector<RRVector> &fOld = this->nodeLoadOldValid ? this->nodeLoadOld : this->nodeLoad;
            for (uint i=0;i<nn;i++)
            {
                RR3Vector fg(fOld.x[i],fOld.y[i],fOld.z[i]);
                RR3Vector fl(fg);
                if (this->localRotations[i].isActive())
                {
                    // f_local = R^T * f_global, see applyLocalRotations().
                    const RRMatrix &R = this->localRotations[i].getR();
                    for (uint c=0;c<3;c++)
                    {
                        fl[c] = R[0][c]*fg[0] + R[1][c]*fg[1] + R[2][c]*fg[2];
                    }
                }
                for (uint c=0;c<3;c++)
                {
                    uint position;
                    if (this->nodeBook.getValue(3*i+c,position))
                    {
                        this->b[position] -= this->hhtAlpha * fl[c];
                    }
                }
            }
        }
    }

    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        RLogger::info("Restoring prestressed nodes\n");
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            this->pModel->setNode(i,nodesBkp[i]);
        }
    }
}

void RSolverStress::solve()
{
    try
    {
        if (this->problemType == R_PROBLEM_STRESS_MODAL)
        {
            this->solveEigenValue();
        }
        else
        {
            this->solveStressStrain();
        }
    }
    catch (const RError &error)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve problem. %s", error.getMessage().toUtf8().constData());
    }
}

void RSolverStress::solveStressStrain()
{
    if (this->computeInitialAcceleration)
    {
        this->solveInitialAcceleration();
    }

    RLogger::info("Solving stress-strain problem.\n");

    try
    {
        RLogger::indent();
        // K + a0*M + a1*C is symmetric positive definite just as K is.
        RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,3);
        RLogger::unindent();
    }
    catch (const RError &error)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve stress-strain problem. %s", error.getMessage().toUtf8().constData());
    }

    this->setDisplacement(this->x);

    if (this->dynamic)
    {
        this->updateVelocityAndAcceleration();
    }
}

bool RSolverStress::findDynamic() const
{
    return (this->problemType == R_PROBLEM_STRESS &&
            this->pModel->getTimeSolver().getEnabled() &&
            this->pModel->getProblemSetup().getStressSetup().getTransientType() == R_STRESS_TRANSIENT_DYNAMIC);
}

void RSolverStress::findTimeIntegrationParameters()
{
    this->timeStepSize = this->pModel->getTimeSolver().getCurrentTimeStepSize();
    if (this->timeStepSize < RConstants::eps)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Dynamic stress analysis requires a positive time step size.");
    }

    this->hhtAlpha = this->pModel->getProblemSetup().getStressSetup().getHhtAlpha();
    this->hhtAlpha = std::min(std::max(this->hhtAlpha,R_STRESS_HHT_ALPHA_MIN_VALUE),R_STRESS_HHT_ALPHA_MAX_VALUE);

    if (this->hhtAlpha < 0.0)
    {
        // Hilber-Hughes-Taylor - second order accurate, unconditionally stable
        // and dissipative in the high frequencies.
        this->newmarkGamma = 0.5 - this->hhtAlpha;
        this->newmarkBeta = std::pow(1.0 - this->hhtAlpha,2) / 4.0;
    }
    else
    {
        // The time-march approximation coefficient is 1.0 for the backward, 0.5
        // for the central and 0.0 for the forward approximation. Newmark with
        // gamma = 1/2 requires beta >= 1/4 for unconditional stability, so the
        // coefficient is clamped to the average-acceleration scheme, the same
        // way as in the acoustic solver.
        this->newmarkGamma = 0.5;
        this->newmarkBeta = std::max(this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient() / 2.0,0.25);
    }
}

void RSolverStress::applyInitialConditions()
{
    uint nn = this->pModel->getNNodes();

    // A new analysis starts from rest unless initial conditions say otherwise.
    this->nodeDisplacement.x.resize(nn,0.0);
    this->nodeDisplacement.y.resize(nn,0.0);
    this->nodeDisplacement.z.resize(nn,0.0);
    this->nodeDisplacement.x.fill(0.0);
    this->nodeDisplacement.y.fill(0.0);
    this->nodeDisplacement.z.fill(0.0);
    this->nodeVelocity.x.resize(nn,0.0);
    this->nodeVelocity.y.resize(nn,0.0);
    this->nodeVelocity.z.resize(nn,0.0);
    this->nodeVelocity.x.fill(0.0);
    this->nodeVelocity.y.fill(0.0);
    this->nodeVelocity.z.fill(0.0);
    this->nodeAcceleration.x.resize(nn,0.0);
    this->nodeAcceleration.y.resize(nn,0.0);
    this->nodeAcceleration.z.resize(nn,0.0);
    this->nodeAcceleration.x.fill(0.0);
    this->nodeAcceleration.y.fill(0.0);
    this->nodeAcceleration.z.fill(0.0);

    RSolverCartesianVector<RRVector> elementValues;
    RSolverCartesianVector<RBVector> setValues;

    this->generateVariableVector(R_VARIABLE_DISPLACEMENT_X,elementValues.x,setValues.x,false,true,false);
    this->generateVariableVector(R_VARIABLE_DISPLACEMENT_Y,elementValues.y,setValues.y,false,true,false);
    this->generateVariableVector(R_VARIABLE_DISPLACEMENT_Z,elementValues.z,setValues.z,false,true,false);
    this->pModel->convertElementToNodeVector(elementValues.x,setValues.x,this->nodeDisplacement.x,true);
    this->pModel->convertElementToNodeVector(elementValues.y,setValues.y,this->nodeDisplacement.y,true);
    this->pModel->convertElementToNodeVector(elementValues.z,setValues.z,this->nodeDisplacement.z,true);

    this->generateVariableVector(R_VARIABLE_STRUCTURAL_VELOCITY_X,elementValues.x,setValues.x,false,true,false);
    this->generateVariableVector(R_VARIABLE_STRUCTURAL_VELOCITY_Y,elementValues.y,setValues.y,false,true,false);
    this->generateVariableVector(R_VARIABLE_STRUCTURAL_VELOCITY_Z,elementValues.z,setValues.z,false,true,false);
    this->pModel->convertElementToNodeVector(elementValues.x,setValues.x,this->nodeVelocity.x,true);
    this->pModel->convertElementToNodeVector(elementValues.y,setValues.y,this->nodeVelocity.y,true);
    this->pModel->convertElementToNodeVector(elementValues.z,setValues.z,this->nodeVelocity.z,true);

    // Constrained directions follow the prescribed displacement and cannot be
    // given a velocity of their own.
    for (uint i=0;i<nn;i++)
    {
        uint nConstrained = this->nodeConstrainedDirections[i];
        if (nConstrained == 0)
        {
            continue;
        }

        RR3Vector vg(this->nodeVelocity.x[i],this->nodeVelocity.y[i],this->nodeVelocity.z[i]);
        RR3Vector vl(vg);

        if (this->localRotations[i].isActive())
        {
            // v_local = R^T * v_global, see applyLocalRotations().
            const RRMatrix &R = this->localRotations[i].getR();
            for (uint c=0;c<3;c++)
            {
                vl[c] = R[0][c]*vg[0] + R[1][c]*vg[1] + R[2][c]*vg[2];
            }
            for (uint c=0;c<nConstrained;c++)
            {
                vl[c] = 0.0;
            }
            for (uint c=0;c<3;c++)
            {
                vg[c] = R[c][0]*vl[0] + R[c][1]*vl[1] + R[c][2]*vl[2];
            }
        }
        else
        {
            for (uint c=0;c<nConstrained;c++)
            {
                vg[c] = 0.0;
            }
        }

        this->nodeVelocity.x[i] = vg[0];
        this->nodeVelocity.y[i] = vg[1];
        this->nodeVelocity.z[i] = vg[2];
    }
}

void RSolverStress::storeOldState()
{
    uint timeStep = this->pModel->getTimeSolver().getCurrentTimeStep();

    // A task group may run the solver several times within one time step. The
    // state at the beginning of the step must then stay as it was.
    if (timeStep == this->oldStateTimeStep && !this->firstRun)
    {
        return;
    }

    uint nn = this->pModel->getNNodes();

    this->nodeLoadOldValid = (this->oldStateTimeStep != RConstants::eod && this->nodeLoad.x.size() == nn);
    if (this->nodeLoadOldValid)
    {
        this->nodeLoadOld = this->nodeLoad;
    }

    this->nodeDisplacementOld = this->nodeDisplacement;
    this->nodeVelocityOld = this->nodeVelocity;
    this->nodeAccelerationOld = this->nodeAcceleration;
    this->oldStateTimeStep = timeStep;

    // At the very beginning the acceleration follows from the equilibrium
    // M*a = f - K*u - C*v. Otherwise a suddenly applied load starts with a
    // wrong acceleration which the non-dissipative scheme never forgets.
    this->computeInitialAcceleration = this->firstRun;
}

void RSolverStress::solveInitialAcceleration()
{
    RLogger::info("Solving initial acceleration.\n");

    uint n = this->M.getNRows();

    // Directions without mass (zero density, massless truss directions) carry
    // no inertia - their initial acceleration is set to zero.
    double maxDiagonal = 0.0;
    RRVector diagonal(n,0.0);
    for (uint i=0;i<n;i++)
    {
        uint position;
        if (this->M.getVector(i).findPosition(i,position))
        {
            diagonal[i] = this->M.getVector(i).getValue(position);
        }
        maxDiagonal = std::max(maxDiagonal,diagonal[i]);
    }

    RRVector a(n,0.0);

    if (maxDiagonal > 0.0)
    {
        for (uint i=0;i<n;i++)
        {
            if (diagonal[i] <= RConstants::eps * maxDiagonal)
            {
                this->M.addValue(i,i,maxDiagonal);
                this->r[i] = 0.0;
            }
        }

        try
        {
            RLogger::indent();
            RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
            matrixSolver.solve(this->M,this->r,a,R_MATRIX_PRECONDITIONER_JACOBI,3);
            RLogger::unindent();
        }
        catch (const RError &error)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve initial acceleration. %s", error.getMessage().toUtf8().constData());
        }

        // The predictor was assembled with zero initial acceleration.
        RRVector Da;
        RSparseMatrix::mlt(this->D,a,Da);
        for (uint i=0;i<n;i++)
        {
            this->b[i] += Da[i];
        }
    }

    // Constrained directions start without acceleration.
    this->setNodeVector(a,this->nodeAccelerationOld,false);
    this->computeInitialAcceleration = false;
}

void RSolverStress::updateVelocityAndAcceleration()
{
    const double gamma = this->newmarkGamma;
    const double beta = this->newmarkBeta;
    const double dt = this->timeStepSize;

    const double a0 = 1.0 / (beta * std::pow(dt,2));
    const double a2 = 1.0 / (beta * dt);
    const double a3 = (1.0 / (2.0 * beta)) - 1.0;
    const double a6 = dt * (1.0 - gamma);
    const double a7 = dt * gamma;

    uint nn = this->pModel->getNNodes();

    this->nodeVelocity.x.resize(nn,0.0);
    this->nodeVelocity.y.resize(nn,0.0);
    this->nodeVelocity.z.resize(nn,0.0);
    this->nodeAcceleration.x.resize(nn,0.0);
    this->nodeAcceleration.y.resize(nn,0.0);
    this->nodeAcceleration.z.resize(nn,0.0);

    RRVector *pU[3] = { &this->nodeDisplacement.x, &this->nodeDisplacement.y, &this->nodeDisplacement.z };
    const RRVector *pUo[3] = { &this->nodeDisplacementOld.x, &this->nodeDisplacementOld.y, &this->nodeDisplacementOld.z };
    const RRVector *pVo[3] = { &this->nodeVelocityOld.x, &this->nodeVelocityOld.y, &this->nodeVelocityOld.z };
    const RRVector *pAo[3] = { &this->nodeAccelerationOld.x, &this->nodeAccelerationOld.y, &this->nodeAccelerationOld.z };
    RRVector *pV[3] = { &this->nodeVelocity.x, &this->nodeVelocity.y, &this->nodeVelocity.z };
    RRVector *pA[3] = { &this->nodeAcceleration.x, &this->nodeAcceleration.y, &this->nodeAcceleration.z };

    for (uint c=0;c<3;c++)
    {
        for (uint i=0;i<nn;i++)
        {
            double uo = (*pUo[c])[i];
            double vo = (*pVo[c])[i];
            double ao = (*pAo[c])[i];
            double u = (*pU[c])[i];

            double a = a0 * (u - uo) - a2 * vo - a3 * ao;
            double v = vo + a6 * ao + a7 * a;

            (*pA[c])[i] = a;
            (*pV[c])[i] = v;
        }
    }
}

void RSolverStress::setNodeVector(const RRVector &v, RSolverCartesianVector<RRVector> &nodeVector, bool prescribed) const
{
    uint nn = this->pModel->getNNodes();

    nodeVector.x.resize(nn,0.0);
    nodeVector.y.resize(nn,0.0);
    nodeVector.z.resize(nn,0.0);

    for (uint i=0;i<nn;i++)
    {
        uint position;
        RR3Vector du(0.0,0.0,0.0);

        for (uint c=0;c<3;c++)
        {
            if (this->nodeBook.getValue(3*i+c,position))
            {
                du[c] = v[position];
            }
            else if (prescribed)
            {
                du[c] = this->nodePrescribedDisplacement[i][c];
            }
        }

        if (this->localRotations[i].isActive())
        {
            this->localRotations[i].rotateResultsVector(du);
        }

        nodeVector.x[i] = du[0];
        nodeVector.y[i] = du[1];
        nodeVector.z[i] = du[2];
    }
}

void RSolverStress::storeNodeVector(RVariableType variableType, const RSolverCartesianVector<RRVector> &nodeVector)
{
    uint nn = this->pModel->getNNodes();

    uint variablePos = this->pModel->findVariable(variableType);
    if (variablePos == RConstants::eod)
    {
        variablePos = this->pModel->addVariable(variableType);

        double vmin = 0.0;
        double vmax = 0.0;
        for (uint i=0;i<nn;i++)
        {
            double v = RR3Vector(nodeVector.x[i],nodeVector.y[i],nodeVector.z[i]).length();
            if (i == 0)
            {
                vmin = vmax = v;
            }
            else
            {
                vmin = std::min(vmin,v);
                vmax = std::max(vmax,v);
            }
        }

        this->pModel->getVariable(variablePos).getVariableData().setMinMaxDisplayValue(vmin,vmax);
    }
    RVariable &variable = this->pModel->getVariable(variablePos);

    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(3,nn);
    for (uint i=0;i<nn;i++)
    {
        variable.setValue(0,i,nodeVector.x[i]);
        variable.setValue(1,i,nodeVector.y[i]);
        variable.setValue(2,i,nodeVector.z[i]);
    }
}

void RSolverStress::lumpMassMatrix(RRMatrix &Me)
{
    // Row sum lumping gives negative masses on corner nodes of quadratic
    // elements. Scaling the diagonal of the consistent matrix to the element
    // mass (HRZ lumping) stays positive for every element type.
    uint nn = Me.getNRows() / 3;

    for (uint d=0;d<3;d++)
    {
        double totalMass = 0.0;
        double diagonalMass = 0.0;
        for (uint m=0;m<nn;m++)
        {
            diagonalMass += Me[3*m+d][3*m+d];
            for (uint n=0;n<nn;n++)
            {
                totalMass += Me[3*m+d][3*n+d];
            }
        }
        double scale = (diagonalMass > 0.0) ? totalMass / diagonalMass : 0.0;
        for (uint m=0;m<nn;m++)
        {
            Me[3*m+d][3*m+d] *= scale;
        }
    }

    for (uint p=0;p<Me.getNRows();p++)
    {
        for (uint q=0;q<Me.getNColumns();q++)
        {
            if (p != q)
            {
                Me[p][q] = 0.0;
            }
        }
    }
}

void RSolverStress::solveEigenValue()
{
    RLogger::info("Solving eigen-value problem.\n");

    REigenValueSolverConf conf;

    if (this->pModel->getProblemSetup().getModalSetup().getMethod() == R_MODAL_MULTIPLE_MODES)
    {
        conf.setMethod(REigenValueSolverConf::SubspaceIteration);
    }
    else
    {
        conf.setMethod(REigenValueSolverConf::InversePowerIteration);
    }
    conf.setNEigenValues(this->pModel->getProblemSetup().getModalSetup().getNModesToExtract());
    conf.setNIterations(this->pModel->getProblemSetup().getModalSetup().getNIterations());
    conf.setSolverCvgValue(this->pModel->getProblemSetup().getModalSetup().getConvergenceValue());
    conf.setOutputFrequency(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getOutputFrequency());
    conf.setOutputFileName(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getOutputFileName());

    REigenValueSolver solver(conf,this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));

    try
    {
        RLogger::indent();
        solver.solve(this->M,this->A,this->d,this->ev);
        RLogger::unindent();
    }
    catch (const RError &error)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve eigen-value problem. %s", error.getMessage().toUtf8().constData());
    }
}

void RSolverStress::setDisplacement(const RRVector &v)
{
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        RR3Vector du(0.0,0.0,0.0);

        for (uint c=0;c<3;c++)
        {
            if (this->nodeBook.getValue(3*i+c,position))
            {
                du[c] = v[position];
            }
            else
            {
                // Constrained direction. The eigen-value problem is homogeneous
                // so prescribed values do not enter the mode shapes.
                du[c] = (this->problemType == R_PROBLEM_STRESS_MODAL) ? 0.0 : this->nodePrescribedDisplacement[i][c];
            }
        }

        if (this->localRotations[i].isActive())
        {
            this->localRotations[i].rotateResultsVector(du);
        }

        if (this->problemType == R_PROBLEM_STRESS_MODAL)
        {
            du[0] += this->nodeInitialDisplacement.x[i];
            du[1] += this->nodeInitialDisplacement.y[i];
            du[2] += this->nodeInitialDisplacement.z[i];
        }

        this->nodeDisplacement.x[i] = du[0];
        this->nodeDisplacement.y[i] = du[1];
        this->nodeDisplacement.z[i] = du[2];
    }
}

void RSolverStress::process()
{

    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        uint modeNum = this->pModel->getProblemSetup().getModalSetup().getMode();

        // The eigen value of K * phi = lambda * M * phi is lambda = omega^2, so
        // the natural frequency is sqrt(lambda) / (2*pi). The modal setup holds
        // a frequency in Hz, not the raw eigen value.
        double eigenValue = this->d[modeNum];
        double frequency = 0.0;

        if (std::isfinite(eigenValue) && eigenValue > 0.0)
        {
            frequency = std::sqrt(eigenValue) / (2.0 * RConstants::pi);
            RLogger::info("Eigen-value = %g, frequency = %g [Hz]\n",eigenValue,frequency);
        }
        else
        {
            RLogger::warning("Mode %u was not resolved - its eigen value is not usable.\n",modeNum+1);
        }

        this->pModel->getProblemSetup().getModalSetup().setFrequency(frequency);

        RRVector v(this->ev.getNColumns(),0.0);
        for (uint j=0;j<this->ev.getNColumns();j++)
        {
            v[j] = this->ev[modeNum][j];
        }
        this->setDisplacement(v);
    }

    // Initialize force vector
    this->nodeForce.x.resize(this->pModel->getNNodes());
    this->nodeForce.y.resize(this->pModel->getNNodes());
    this->nodeForce.z.resize(this->pModel->getNNodes());

    this->nodeForce.x.fill(0.0);
    this->nodeForce.y.fill(0.0);
    this->nodeForce.z.fill(0.0);

    // Initialize stress vectors
    for (uint i=0;i<6;i++)
    {
        this->elementStress[i].resize(this->pModel->getNElements(),0.0);
        this->elementStress[i].fill(0.0);
    }
    this->elementNormalStress.resize(this->pModel->getNElements(),0.0);
    this->elementShearStress.resize(this->pModel->getNElements(),0.0);
    this->elementVonMisses.resize(this->pModel->getNElements(),0.0);

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);
        double lineCrossArea = line.getCrossArea();

        if (lineCrossArea == 0.0)
        {
            continue;
        }

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = line.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                double QeN = 0.0;

                RRMatrix Be(3*element.size(),1);
                RRMatrix BeT(1,3*element.size());

                double E = this->elementElasticityModulus[elementID];
                // Axial stiffness of the truss - used for the nodal force only.
                double De = E * lineCrossArea;

                RRMatrix Rl;
                RRVector tl;
                element.findTransformationMatrix(this->pModel->getNodes(),Rl,tl);
                Rl.invert();

                RRVector lxe(element.size(),0.0);
                for (uint k=0;k<element.size();k++)
                {
                    RR3Vector xg(this->nodeDisplacement.x[element.getNodeId(k)],
                                 this->nodeDisplacement.y[element.getNodeId(k)],
                                 this->nodeDisplacement.z[element.getNodeId(k)]);
                    RR3Vector xl;
                    RRMatrix::mlt(Rl,xg,xl);
                    lxe[k] = xl[0];
                }

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];
                RRVector feThermal(element.size()*3,0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    RtT.transpose(Rt);

                    Be.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        Be[3*m+0][0] = Rt[3*m+0][0]*dN[m][0]*J[0][0];
                        Be[3*m+1][0] = Rt[3*m+1][0]*dN[m][0]*J[0][0];
                        Be[3*m+2][0] = Rt[3*m+2][0]*dN[m][0]*J[0][0];
                    }
                    BeT.transpose(Be);

                    // Same stiffness as the one assembled in prepare().
                    RRMatrix BeScaled(Be);
                    BeScaled *= De * detJ * shapeFunc.getW();
                    RRMatrix::mlt(BeScaled,BeT,Ke,true);

                    // Same thermal load as the one assembled in prepare().
                    double thermalFactor = De
                                         * this->elementThermalExpansion[elementID]
                                         * dT
                                         * detJ
                                         * shapeFunc.getW();
                    for (uint m=0;m<3*element.size();m++)
                    {
                        feThermal[m] += thermalFactor * Be[m][0];
                    }

                    double integValue = 1.0/double(nInp);

                    // Element level stress. The axial strain follows from the
                    // local axial displacements; the stress is E*(eps - alpha*dT)
                    // and must not carry the cross area, which belongs to the
                    // stiffness only.
                    double axialStrain = 0.0;
                    for (uint m=0;m<element.size();m++)
                    {
                        axialStrain += dN[m][0]*J[0][0] * lxe[m];
                    }
                    QeN += E * (axialStrain - this->elementThermalExpansion[elementID] * dT) * integValue;
                }

                // Nodal force is the internal elastic force K*u - f_thermal.
                // There is no acceleration state in the formulation, so no
                // inertia term.
                RRMatrix::mlt(Ke,xe,fe);
                for (uint n=0;n<fe.size();n++)
                {
                    fe[n] -= feThermal[n];
                }

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // A truss carries an axial stress only, along its local x axis.
                this->elementStress[0][elementID] = QeN;
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = 0.0;
                this->elementVonMisses[elementID] = std::fabs(QeN);
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);
        double surfaceThickness = surface.getThickness();

        if (surfaceThickness == 0.0)
        {
            continue;
        }

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = surface.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                RRVector Qe(3,0.0);
                double QeN = 0.0;
                double QeS = 0.0;
                double QeVM = 0.0;

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*2,3);
                RRMatrix BeT(3,element.size()*2);
                RRMatrix BeD(element.size()*2,3);
                RRMatrix Ket(element.size()*2,element.size()*2);
                RRMatrix KeRt(element.size()*2,element.size()*2);
                RRVector fet(element.size()*2);

                RRMatrix De(3,3,0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1-v;   De[0][1] = v;
                De[1][0] = v;     De[1][1] = 1-v;
                De[2][2] = (1-2*v)/2;
                De *= E/((1+v)*(1-2*v));

                RRMatrix Rl;
                RRVector tl;
                element.findTransformationMatrix(this->pModel->getNodes(),Rl,tl);
                Rl.invert();

                RRVector lxe(element.size()*2);
                for (uint k=0;k<element.size();k++)
                {
                    RR3Vector xg(this->nodeDisplacement.x[element.getNodeId(k)],
                                 this->nodeDisplacement.y[element.getNodeId(k)],
                                 this->nodeDisplacement.z[element.getNodeId(k)]);
                    RR3Vector xl;
                    RRMatrix::mlt(Rl,xg,xl);
                    lxe[2*k+0] = xl[0];
                    lxe[2*k+1] = xl[1];
                }

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];
                // Plane strain - see the thermal load in prepare().
                double thermalStrain = (1.0 + v) * this->elementThermalExpansion[elementID] * dT;
                // Thermal load in the global frame.
                RRVector feThermal(element.size()*3,0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    Rt = findSurfaceRotation(Rt);
                    RtT.transpose(Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        Be[2*m][0] = B[m][0];   Be[2*m+1][0] = 0.0;
                        Be[2*m][1] = 0.0;       Be[2*m+1][1] = B[m][1];
                        Be[2*m][2] = B[m][1];   Be[2*m+1][2] = B[m][0];
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    // Same stiffness as assembled in prepare(), accumulated over
                    // the integration points.
                    Ket *= surfaceThickness * detJ * shapeFunc.getW();
                    RRMatrix::mlt(Rt,Ket,KeRt);
                    RRMatrix::mlt(KeRt,RtT,Ke,true);

                    double integValue = 1.0/double(nInp);

                    // Element level stress D*B*u.
                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<3;n++)
                        {
                            Qe[n] += BeD[2*m+0][n] * lxe[2*m+0] * integValue
                                  +  BeD[2*m+1][n] * lxe[2*m+1] * integValue;
                        }

                        double fl0 = 0.0;
                        double fl1 = 0.0;
                        for (uint n=0;n<2;n++)
                        {
                            fl0 += thermalStrain * BeD[2*m+0][n] * surfaceThickness * detJ * shapeFunc.getW();
                            fl1 += thermalStrain * BeD[2*m+1][n] * surfaceThickness * detJ * shapeFunc.getW();
                        }
                        for (uint c=0;c<3;c++)
                        {
                            feThermal[3*m+c] += Rt[3*m+c][2*m+0]*fl0 + Rt[3*m+c][2*m+1]*fl1;
                        }
                    }
                }

                // Thermal stress -D*eps_th, eps_th = thermalStrain * (1,1,0).
                // The stress is uniform over the element, so it is subtracted
                // once and not per node - the shape function derivatives of
                // all nodes sum up to zero.
                for (uint n=0;n<3;n++)
                {
                    Qe[n] -= thermalStrain * (De[n][0] + De[n][1]);
                }

                // Nodal force is the internal elastic force K*u - f_thermal.
                // There is no acceleration state in the formulation, so no
                // inertia term.
                RRMatrix::mlt(Ke,xe,fe);
                for (uint n=0;n<fe.size();n++)
                {
                    fe[n] -= feThermal[n];
                }

                QeN = std::sqrt(Qe[0] * Qe[0] + Qe[1] * Qe[1] - Qe[0] * Qe[1]);
                QeS = std::sqrt(3.0 * Qe[2] * Qe[2]);
                // Von Mises combines the normal and the shear invariant in
                // quadrature, not by adding them.
                QeVM = std::sqrt(QeN * QeN + QeS * QeS);

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // In-plane components, in the local element frame.
                this->elementStress[0][elementID] = Qe[0];
                this->elementStress[1][elementID] = Qe[1];
                this->elementStress[5][elementID] = Qe[2];
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = QeS;
                this->elementVonMisses[elementID] = QeVM;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }

    // Process volume elements.
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
                uint elementID = volume.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                RRVector Qe(6,0.0);

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*3,6);
                RRMatrix BeT(6,element.size()*3);
                RRMatrix BeD(element.size()*3,6);
                RRMatrix Ket(element.size()*3,element.size()*3);

                RRMatrix De(6,6,0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1-v;   De[0][1] = v;     De[0][2] = v;
                De[1][0] = v;     De[1][1] = 1-v;   De[1][2] = v;
                De[2][0] = v;     De[2][1] = v;     De[2][2] = 1-v;
                De[3][3] = De[4][4] = De[5][5] = (1-2*v)/2;
                De *= E/((1+v)*(1-2*v));

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];
                double thermalStrain = this->elementThermalExpansion[elementID] * dT;
                RRVector feThermal(element.size()*3,0.0);

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        Be[3*m+0][0] = B[m][0];   Be[3*m+1][0] = 0.0;       Be[3*m+2][0] = 0.0;
                        Be[3*m+0][1] = 0.0;       Be[3*m+1][1] = B[m][1];   Be[3*m+2][1] = 0.0;
                        Be[3*m+0][2] = 0.0;       Be[3*m+1][2] = 0.0;       Be[3*m+2][2] = B[m][2];
                        Be[3*m+0][3] = 0.0;       Be[3*m+1][3] = B[m][2];   Be[3*m+2][3] = B[m][1];
                        Be[3*m+0][4] = B[m][2];   Be[3*m+1][4] = 0.0;       Be[3*m+2][4] = B[m][0];
                        Be[3*m+0][5] = B[m][1];   Be[3*m+1][5] = B[m][0];   Be[3*m+2][5] = 0.0;
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    for (uint m=0;m<3*element.size();m++)
                    {
                        for (uint n=0;n<3*element.size();n++)
                        {
                            // Stiffness matrix
                            Ke[m][n] += Ket[m][n] * detJ * shapeFunc.getW();
                        }
                    }

                    double integValue = 1.0/double(nInp);

                    // Element level stress.
                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<6;n++)
                        {
                            Qe[n] += BeD[3*m+0][n] * this->nodeDisplacement.x[element.getNodeId(m)] * integValue
                                  +  BeD[3*m+1][n] * this->nodeDisplacement.y[element.getNodeId(m)] * integValue
                                  +  BeD[3*m+2][n] * this->nodeDisplacement.z[element.getNodeId(m)] * integValue;
                        }
                        // Thermal load B^T * D * eps_th, as in prepare().
                        for (uint n=0;n<3;n++)
                        {
                            feThermal[3*m+0] += thermalStrain * BeD[3*m+0][n] * detJ * shapeFunc.getW();
                            feThermal[3*m+1] += thermalStrain * BeD[3*m+1][n] * detJ * shapeFunc.getW();
                            feThermal[3*m+2] += thermalStrain * BeD[3*m+2][n] * detJ * shapeFunc.getW();
                        }
                    }
                }

                // Thermal stress -D*eps_th, eps_th = thermalStrain * (1,1,1,0,0,0).
                // The stress is uniform over the element, so it is subtracted
                // once and not per node - the shape function derivatives of
                // all nodes sum up to zero.
                for (uint n=0;n<6;n++)
                {
                    Qe[n] -= thermalStrain * (De[n][0] + De[n][1] + De[n][2]);
                }

                // Nodal force is the internal elastic force K*u - f_thermal.
                // There is no acceleration state in the formulation, so no
                // inertia term.
                RRMatrix::mlt(Ke,xe,fe);
                for (uint n=0;n<fe.size();n++)
                {
                    fe[n] -= feThermal[n];
                }

                double QeN = std::sqrt(Qe[0]*Qe[0] + Qe[1]*Qe[1] + Qe[2]*Qe[2] - (Qe[0]*Qe[1] + Qe[1]*Qe[2] + Qe[2]*Qe[0]));
                double QeS = std::sqrt(3.0 * (Qe[3]*Qe[3] + Qe[4]*Qe[4] + Qe[5]*Qe[5]));
                // Von Mises combines the normal and the shear invariant in
                // quadrature, not by adding them.
                double QeVM = std::sqrt(QeN*QeN + QeS*QeS);

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // Stress components in global coordinates.
                for (uint n=0;n<6;n++)
                {
                    this->elementStress[n][elementID] = Qe[n];
                }
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = QeS;
                this->elementVonMisses[elementID] = QeVM;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }
}

void RSolverStress::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Displacement
    uint displacementPos = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
    if (displacementPos == RConstants::eod)
    {
        displacementPos = this->pModel->addVariable(R_VARIABLE_DISPLACEMENT);

        double umin = 0.0;
        double umax = 0.0;
        for (uint i=0;i<this->nodeDisplacement.x.size();i++)
        {
            double u = RR3Vector(this->nodeDisplacement.x[i],
                                 this->nodeDisplacement.y[i],
                                 this->nodeDisplacement.z[i]).length();
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

        this->pModel->getVariable(displacementPos).getVariableData().setMinMaxDisplayValue(umin,umax);
    }
    RVariable &displacement =  this->pModel->getVariable(displacementPos);

    displacement.setApplyType(R_VARIABLE_APPLY_NODE);
    displacement.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        displacement.setValue(0,i,this->nodeDisplacement.x[i]);
        displacement.setValue(1,i,this->nodeDisplacement.y[i]);
        displacement.setValue(2,i,this->nodeDisplacement.z[i]);
    }

    if (this->dynamic)
    {
        // Velocity and acceleration - also carry the Newmark state over to the
        // next time step.
        this->storeNodeVector(R_VARIABLE_STRUCTURAL_VELOCITY,this->nodeVelocity);
        this->storeNodeVector(R_VARIABLE_STRUCTURAL_ACCELERATION,this->nodeAcceleration);
    }

    // VonMises Stress
    uint vonMisesStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_VON_MISES);
    if (vonMisesStressPos == RConstants::eod)
    {
        vonMisesStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_VON_MISES);

        this->pModel->getVariable(vonMisesStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementVonMisses),
                    RStatistics::findMaximumValue(this->elementVonMisses));
    }
    RVariable &vonMisesStress =  this->pModel->getVariable(vonMisesStressPos);

    vonMisesStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    vonMisesStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        vonMisesStress.setValue(0,i,this->elementVonMisses[i]);
    }

    // Stress components. Volume elements report them in global coordinates,
    // surface and line elements in their own local element frame.
    const RVariableType stressComponentTypes[6] =
    {
        R_VARIABLE_STRESS_X,
        R_VARIABLE_STRESS_Y,
        R_VARIABLE_STRESS_Z,
        R_VARIABLE_STRESS_YZ,
        R_VARIABLE_STRESS_XZ,
        R_VARIABLE_STRESS_XY
    };

    for (uint c=0;c<6;c++)
    {
        uint stressComponentPos = this->pModel->findVariable(stressComponentTypes[c]);
        if (stressComponentPos == RConstants::eod)
        {
            stressComponentPos = this->pModel->addVariable(stressComponentTypes[c]);

            this->pModel->getVariable(stressComponentPos).getVariableData().setMinMaxDisplayValue(
                        RStatistics::findMinimumValue(this->elementStress[c]),
                        RStatistics::findMaximumValue(this->elementStress[c]));
        }
        RVariable &stressComponent = this->pModel->getVariable(stressComponentPos);

        stressComponent.setApplyType(R_VARIABLE_APPLY_ELEMENT);
        stressComponent.resize(1,this->pModel->getNElements());
        for (uint i=0;i<this->pModel->getNElements();i++)
        {
            stressComponent.setValue(0,i,this->elementStress[c][i]);
        }
    }

    // Normal Stress
    uint normalStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_NORMAL);
    if (normalStressPos == RConstants::eod)
    {
        normalStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_NORMAL);

        this->pModel->getVariable(normalStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementNormalStress),
                    RStatistics::findMaximumValue(this->elementNormalStress));
    }
    RVariable &normalStress =  this->pModel->getVariable(normalStressPos);

    normalStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    normalStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        normalStress.setValue(0,i,this->elementNormalStress[i]);
    }

    // Shear Stress
    uint shearStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_SHEAR);
    if (shearStressPos == RConstants::eod)
    {
        shearStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_SHEAR);

        this->pModel->getVariable(shearStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementShearStress),
                    RStatistics::findMaximumValue(this->elementShearStress));
    }
    RVariable &shearStress =  this->pModel->getVariable(shearStressPos);

    shearStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    shearStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        shearStress.setValue(0,i,this->elementShearStress[i]);
    }

    // Force
    uint forcePos = this->pModel->findVariable(R_VARIABLE_FORCE);
    if (forcePos == RConstants::eod)
    {
        forcePos = this->pModel->addVariable(R_VARIABLE_FORCE);

        double fmin = 0.0;
        double fmax = 0.0;
        for (uint i=0;i<this->nodeForce.x.size();i++)
        {
            double f = RR3Vector(this->nodeForce.x[i],
                                 this->nodeForce.y[i],
                                 this->nodeForce.z[i]).length();
            if (i == 0)
            {
                fmin = fmax = f;
            }
            else
            {
                fmin = std::min(fmin,f);
                fmax = std::max(fmax,f);
            }
        }

        this->pModel->getVariable(forcePos).getVariableData().setMinMaxDisplayValue(fmin,fmax);
    }
    RVariable &force =  this->pModel->getVariable(forcePos);

    force.setApplyType(R_VARIABLE_APPLY_NODE);
    force.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        force.setValue(0,i,this->nodeForce.x[i]);
        force.setValue(1,i,this->nodeForce.y[i]);
        force.setValue(2,i,this->nodeForce.z[i]);
    }

    RLogger::unindent();
}

void RSolverStress::statistics()
{
    this->printStats(R_VARIABLE_DISPLACEMENT);
    if (this->dynamic)
    {
        this->printStats(R_VARIABLE_STRUCTURAL_VELOCITY);
        this->printStats(R_VARIABLE_STRUCTURAL_ACCELERATION);
    }
    this->printStats(R_VARIABLE_STRESS_VON_MISES);
    this->printStats(R_VARIABLE_FORCE);
    this->processMonitoringPoints();
}

bool RSolverStress::isComponentEnabled(const RBoundaryCondition &bc, RVariableType variableType)
{
    uint componentPosition = bc.findComponentPosition(variableType);
    if (componentPosition == RConstants::eod)
    {
        return false;
    }
    return bc.getComponent(componentPosition).getEnabled();
}

double RSolverStress::findComponentValue(const RBoundaryCondition &bc, RVariableType variableType) const
{
    uint componentPosition = bc.findComponentPosition(variableType);
    if (componentPosition == RConstants::eod)
    {
        return 0.0;
    }
    return bc.getComponent(componentPosition).get(this->pModel->getTimeSolver().getCurrentTime());
}

void RSolverStress::updateLocalRotations()
{
    // The local frame of a node follows from the constraints acting on it, so
    // it is built in generateLocalConstraints() instead of from the geometry of
    // a single boundary condition.
}

namespace
{

//! One displacement constraint acting on a node: direction . u = value, with
//! the direction given in global coordinates and of unit length.
struct RNodeConstraint
{
    RR3Vector direction;
    double value;
};

//! Return the boundary condition of an entity which carries a local direction,
//! or null when there is none.
const RBoundaryCondition *findLocalDirectionBc(const RElementGroup &rElementGroup, RProblemType problemType)
{
    for (uint i=0;i<rElementGroup.getNBoundaryConditions();i++)
    {
        const RBoundaryCondition &bc = rElementGroup.getBoundaryCondition(i);
        if ((RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType) && bc.getHasLocalDirection())
        {
            return &bc;
        }
    }
    return nullptr;
}

}

void RSolverStress::generateLocalConstraints()
{
    uint nNodes = this->pModel->getNNodes();

    this->localRotations.resize(nNodes);
    this->nodeConstrainedDirections.resize(nNodes,0);
    this->nodeConstrainedDirections.fill(0);
    this->nodePrescribedDisplacement.resize(nNodes,3,0.0);
    this->nodePrescribedDisplacement.fill(0.0);

    // Directions taken from the geometry. A node lying on several entities
    // collects the average of what they give it, as before.
    std::vector<RR3Vector> nodeNormal(nNodes,RR3Vector(0.0,0.0,0.0));
    RBVector nodeNormalSet(nNodes,false);
    std::vector<RR3Vector> nodeLineDirection(nNodes,RR3Vector(0.0,0.0,0.0));
    RBVector nodeLineDirectionSet(nNodes,false);

    // Surfaces - the averaged element normals, or the entered direction.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rSurface,this->problemType);
        if (!pBc)
        {
            continue;
        }

        bool useEntered = pBc->getExplicitLocalDirection();
        RR3Vector entered = pBc->getLocalDirection();

        if (useEntered && entered.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of surface entity \'%s\' has zero length.",
                         rSurface.getName().toUtf8().constData());
        }

        for (uint j=0;j<rSurface.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rSurface.get(j));
            RR3Vector direction;

            if (useEntered)
            {
                direction = entered;
            }
            else if (!rElement.findNormal(this->pModel->getNodes(),direction[0],direction[1],direction[2]))
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Could not calculate element normal for element# = %u.",rSurface.get(j));
            }

            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                nodeNormal[nodeId][0] += direction[0];
                nodeNormal[nodeId][1] += direction[1];
                nodeNormal[nodeId][2] += direction[2];
                nodeNormalSet[nodeId] = true;
            }
        }
    }

    // Points - there is no geometry to follow, the entered direction is used.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        const RPoint &rPoint = this->pModel->getPoint(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rPoint,this->problemType);
        if (!pBc)
        {
            continue;
        }

        RR3Vector direction = pBc->getLocalDirection();
        if (direction.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of point entity \'%s\' has zero length.",
                         rPoint.getName().toUtf8().constData());
        }

        for (uint j=0;j<rPoint.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rPoint.get(j));
            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                nodeNormal[nodeId] = direction;
                nodeNormalSet[nodeId] = true;
            }
        }
    }

    // Lines - the element direction, or the entered direction.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        const RLine &rLine = this->pModel->getLine(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rLine,this->problemType);
        if (!pBc)
        {
            continue;
        }

        bool useEntered = pBc->getExplicitLocalDirection();
        RR3Vector entered = pBc->getLocalDirection();

        if (useEntered && entered.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of line entity \'%s\' has zero length.",
                         rLine.getName().toUtf8().constData());
        }

        for (uint j=0;j<rLine.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rLine.get(j));
            if (rElement.size() < 2)
            {
                continue;
            }

            RR3Vector direction;
            if (useEntered)
            {
                direction = entered;
            }
            else
            {
                RR3Vector::subtract(this->pModel->getNode(rElement.getNodeId(1)).toVector(),
                                    this->pModel->getNode(rElement.getNodeId(0)).toVector(),
                                    direction);
            }
            if (direction.normalize() < RConstants::eps)
            {
                continue;
            }

            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                // Keep the accumulated direction consistently oriented, so that
                // a polyline does not average itself away.
                double sign = 1.0;
                if (nodeLineDirectionSet[nodeId] && RR3Vector::dot(nodeLineDirection[nodeId],direction) < 0.0)
                {
                    sign = -1.0;
                }
                nodeLineDirection[nodeId][0] += sign*direction[0];
                nodeLineDirection[nodeId][1] += sign*direction[1];
                nodeLineDirection[nodeId][2] += sign*direction[2];
                nodeLineDirectionSet[nodeId] = true;
            }
        }
    }

    for (uint i=0;i<nNodes;i++)
    {
        if (nodeNormalSet[i])
        {
            nodeNormalSet[i] = (nodeNormal[i].normalize() > RConstants::eps);
        }
        if (nodeLineDirectionSet[i])
        {
            nodeLineDirectionSet[i] = (nodeLineDirection[i].normalize() > RConstants::eps);
        }
    }

    // Collect the constraints of every node, in global coordinates.
    std::vector<std::vector<RNodeConstraint> > nodeConstraints(nNodes);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        REntityGroupType entityType = this->pModel->getEntityGroupType(i);
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (!(RBoundaryCondition::getProblemTypeMask(bc.getType()) & R_PROBLEM_STRESS))
            {
                continue;
            }
            if (bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT &&
                bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT_NORMAL &&
                bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER)
            {
                continue;
            }

            for (uint k=0;k<pElementGroup->size();k++)
            {
                const RElement &rElement = this->pModel->getElement(pElementGroup->get(k));

                for (uint l=0;l<rElement.size();l++)
                {
                    uint nodeId = rElement.getNodeId(l);
                    std::vector<RNodeConstraint> &constraints = nodeConstraints[nodeId];

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT)
                    {
                        // Each switched on component holds one global direction.
                        const RVariableType components[3] =
                        {
                            R_VARIABLE_DISPLACEMENT_X,
                            R_VARIABLE_DISPLACEMENT_Y,
                            R_VARIABLE_DISPLACEMENT_Z
                        };
                        for (uint c=0;c<3;c++)
                        {
                            if (!RSolverStress::isComponentEnabled(bc,components[c]))
                            {
                                continue;
                            }
                            RNodeConstraint constraint;
                            constraint.direction = RR3Vector(0.0,0.0,0.0);
                            constraint.direction[c] = 1.0;
                            constraint.value = this->findComponentValue(bc,components[c]);
                            constraints.push_back(constraint);
                        }
                        continue;
                    }

                    double value = this->findComponentValue(bc,R_VARIABLE_DISPLACEMENT);

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER &&
                        entityType == R_ENTITY_GROUP_LINE)
                    {
                        // Free to slide along the line, held across it. The
                        // prescribed value has no single direction to act in
                        // here, so both held directions are held at zero.
                        if (!nodeLineDirectionSet[nodeId])
                        {
                            continue;
                        }
                        RRMatrix frame;
                        nodeLineDirection[nodeId].findRotationMatrix(frame);

                        for (uint c=1;c<3;c++)
                        {
                            RNodeConstraint constraint;
                            constraint.direction = RR3Vector(frame[0][c],frame[1][c],frame[2][c]);
                            constraint.value = 0.0;
                            constraints.push_back(constraint);
                        }
                        continue;
                    }

                    if (!nodeNormalSet[nodeId])
                    {
                        continue;
                    }

                    RNodeConstraint constraint;
                    constraint.direction = nodeNormal[nodeId];
                    constraint.value = value;
                    constraints.push_back(constraint);

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT_NORMAL)
                    {
                        // Held in the two directions across the normal as well.
                        RRMatrix frame;
                        nodeNormal[nodeId].findRotationMatrix(frame);

                        for (uint c=1;c<3;c++)
                        {
                            RNodeConstraint tangential;
                            tangential.direction = RR3Vector(frame[0][c],frame[1][c],frame[2][c]);
                            tangential.value = 0.0;
                            constraints.push_back(tangential);
                        }
                    }
                }
            }
        }
    }

    // Reduce the constraints of each node to an orthonormal set and build the
    // node frame from it.
    for (uint i=0;i<nNodes;i++)
    {
        const std::vector<RNodeConstraint> &constraints = nodeConstraints[i];
        if (constraints.empty())
        {
            this->localRotations[i].deactivate();
            continue;
        }

        RR3Vector axis[3];
        double prescribed[3] = { 0.0, 0.0, 0.0 };
        uint nConstrained = 0;

        double valueScale = 0.0;
        for (uint j=0;j<constraints.size();j++)
        {
            valueScale = std::max(valueScale,std::fabs(constraints[j].value));
        }

        // Every record is examined, including those beyond the third one - a
        // direction which is already held may still disagree about its value.
        for (uint j=0;j<constraints.size();j++)
        {
            RR3Vector w(constraints[j].direction);
            double wValue = constraints[j].value;

            // Gram-Schmidt against what has been accepted, carrying the
            // prescribed values through the same operations.
            for (uint l=0;l<nConstrained;l++)
            {
                double overlap = RR3Vector::dot(constraints[j].direction,axis[l]);
                w[0] -= overlap*axis[l][0];
                w[1] -= overlap*axis[l][1];
                w[2] -= overlap*axis[l][2];
                wValue -= overlap*prescribed[l];
            }

            double norm = w.length();
            if (norm > 1.0e-8 && nConstrained < 3)
            {
                axis[nConstrained] = RR3Vector(w[0]/norm,w[1]/norm,w[2]/norm);
                prescribed[nConstrained] = wValue/norm;
                nConstrained++;
            }
            else if (std::fabs(wValue) > 1.0e-8*std::max(valueScale,1.0e-12))
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Node %u carries displacement constraints which contradict each other. "
                             "Check that the entities meeting at this node do not prescribe different "
                             "displacements in the same direction.",i);
            }
        }

        // Everything accepted so far is held; what follows only completes the
        // frame and stays free.
        uint nHeld = nConstrained;

        for (uint c=0;c<3 && nConstrained<3;c++)
        {
            RR3Vector candidate(0.0,0.0,0.0);
            candidate[c] = 1.0;

            for (uint l=0;l<nConstrained;l++)
            {
                double overlap = RR3Vector::dot(candidate,axis[l]);
                candidate[0] -= overlap*axis[l][0];
                candidate[1] -= overlap*axis[l][1];
                candidate[2] -= overlap*axis[l][2];
            }

            double norm = candidate.length();
            if (norm > 1.0e-8)
            {
                axis[nConstrained] = RR3Vector(candidate[0]/norm,candidate[1]/norm,candidate[2]/norm);
                prescribed[nConstrained] = 0.0;
                nConstrained++;
            }
        }

        this->nodeConstrainedDirections[i] = nHeld;

        RRMatrix R(3,3,0.0);
        for (uint r=0;r<3;r++)
        {
            for (uint c=0;c<3;c++)
            {
                R[r][c] = axis[c][r];
            }
        }

        for (uint c=0;c<3;c++)
        {
            this->nodePrescribedDisplacement[i][c] = prescribed[c];
        }

        // A frame which is already the global one needs no rotation at all.
        bool isIdentity = true;
        for (uint r=0;r<3 && isIdentity;r++)
        {
            for (uint c=0;c<3;c++)
            {
                double expected = (r == c) ? 1.0 : 0.0;
                if (std::fabs(R[r][c]-expected) > 1.0e-12)
                {
                    isIdentity = false;
                    break;
                }
            }
        }

        if (isIdentity)
        {
            this->localRotations[i].deactivate();
        }
        else
        {
            this->localRotations[i].activate(R);
        }
    }
}

void RSolverStress::generateNodeBook()
{
    this->nodeBook.resize(this->pModel->getNNodes()*3);
    this->nodeBook.initialize();

    // The held directions of a node are the leading directions of its frame.
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        for (uint j=0;j<this->nodeConstrainedDirections[i];j++)
        {
            this->nodeBook.disable(3*i+j,true);
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
            this->nodeBook.disable(3*i+0,true);
            this->nodeBook.disable(3*i+1,true);
            this->nodeBook.disable(3*i+2,true);
        }
    }
}

void RSolverStress::assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, AssemblyBuffer &buffer)
{
    const RElement &rElement = this->pModel->getElement(elementID);

    uint ne = 3*rElement.size();

    RRMatrix Ae(ne,ne);
    RRMatrix Be(ne,ne);
    RRMatrix De(ne,ne);
    RRVector be(ne);
    RRVector re(ne);

    Ae.fill(0.0);
    Be.fill(0.0);
    De.fill(0.0);
    be.fill(0.0);
    re.fill(0.0);

    if (this->dynamic)
    {
        // Newmark time integration of M*u'' + C*u' + K*u = f with the
        // Hilber-Hughes-Taylor modification. Everything is assembled in the
        // global frame and rotated to the node frames afterwards.
        //
        //   (a0*M + (1+alpha)*(K + a1*C)) * u(n+1) =
        //       (1+alpha)*f(n+1) - alpha*f(n)
        //     + M*(a0*u + a2*v + a3*a)
        //     + (1+alpha)*C*(a1*u + a4*v + a5*a)
        //     + alpha*(C*v + K*u)
        //
        // with u, v, a taken at the beginning of the step. The -alpha*f(n) term
        // is added to the assembled system in prepare().
        const double gamma = this->newmarkGamma;
        const double beta = this->newmarkBeta;
        const double alpha = this->hhtAlpha;
        const double dt = this->timeStepSize;

        const double a0 = 1.0 / (beta * std::pow(dt,2));
        const double a1 = gamma / (beta * dt);
        const double a2 = 1.0 / (beta * dt);
        const double a3 = (1.0 / (2.0 * beta)) - 1.0;
        const double a4 = (gamma / beta) - 1.0;
        const double a5 = (1.0 / 2.0) * dt * ((gamma / beta) - 2.0);

        RRMatrix Mm(Me);
        if (this->pModel->getProblemSetup().getStressSetup().getLumpedMass())
        {
            RSolverStress::lumpMassMatrix(Mm);
        }

        const double aM = this->elementRayleighMass[elementID];
        const double aK = this->elementRayleighStiffness[elementID];

        RRVector ue(ne), ve(ne), ae(ne);
        for (uint n=0;n<rElement.size();n++)
        {
            uint nodeID = rElement.getNodeId(n);
            ue[3*n+0] = this->nodeDisplacementOld.x[nodeID];
            ue[3*n+1] = this->nodeDisplacementOld.y[nodeID];
            ue[3*n+2] = this->nodeDisplacementOld.z[nodeID];
            ve[3*n+0] = this->nodeVelocityOld.x[nodeID];
            ve[3*n+1] = this->nodeVelocityOld.y[nodeID];
            ve[3*n+2] = this->nodeVelocityOld.z[nodeID];
            ae[3*n+0] = this->nodeAccelerationOld.x[nodeID];
            ae[3*n+1] = this->nodeAccelerationOld.y[nodeID];
            ae[3*n+2] = this->nodeAccelerationOld.z[nodeID];
        }

        for (uint p=0;p<ne;p++)
        {
            be[p] = (1.0 + alpha) * fe[p];
            re[p] = fe[p];
            for (uint q=0;q<ne;q++)
            {
                double Mpq = Mm[p][q];
                double Kpq = Ke[p][q];
                double Cpq = aM * Mpq + aK * Kpq;

                Ae[p][q] = a0 * Mpq + (1.0 + alpha) * (Kpq + a1 * Cpq);

                be[p] += Mpq * (a0 * ue[q] + a2 * ve[q] + a3 * ae[q])
                      +  (1.0 + alpha) * Cpq * (a1 * ue[q] + a4 * ve[q] + a5 * ae[q])
                      +  alpha * (Cpq * ve[q] + Kpq * ue[q]);

                if (this->computeInitialAcceleration)
                {
                    Be[p][q] = Mpq;
                    De[p][q] = a3 * Mpq + (1.0 + alpha) * a5 * Cpq;
                    re[p] -= Kpq * ue[q] + Cpq * ve[q];
                }
            }
        }

        for (uint m=0;m<rElement.size();m++)
        {
            uint nodeID = rElement.getNodeId(m);
            buffer.f[3*nodeID+0] += fe[3*m+0];
            buffer.f[3*nodeID+1] += fe[3*m+1];
            buffer.f[3*nodeID+2] += fe[3*m+2];
        }
    }
    else if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        Ae = Ke;
        Be = Me;
        be = fe;
    }
    else
    {
        // Static analysis, or a quasi-static sequence of static solutions.
        Ae = Ke;
        be = fe;
    }
    this->applyLocalRotations(elementID,Ae);
    this->applyLocalRotations(elementID,Be);
    this->applyLocalRotations(elementID,be);
    if (this->computeInitialAcceleration)
    {
        this->applyLocalRotations(elementID,De);
        this->applyLocalRotations(elementID,re);
    }

    // Apply explicit boundary conditions. Ae and be are in the frame of each
    // node by now, so the prescribed values have to be taken in that same
    // frame - which is what nodePrescribedDisplacement holds. The eigen-value
    // problem is homogeneous, so nothing is prescribed there.
    if (this->problemType != R_PROBLEM_STRESS_MODAL)
    {
        for (uint m=0;m<rElement.size();m++)
        {
            uint position;
            uint nodeID = rElement.getNodeId(m);
            for (uint c=0;c<3;c++)
            {
                if (this->nodeBook.getValue(3*nodeID+c,position))
                {
                    continue;
                }
                double u = this->nodePrescribedDisplacement[nodeID][c];
                if (std::fabs(u) <= 0.0)
                {
                    continue;
                }
                for (uint n=0;n<ne;n++)
                {
                    be[n] -= Ae[n][3*m+c] * u;
                }
            }
        }
    }

    // Assembly final matrix system
    uint dims = 3;
    const bool assembleMass = (this->problemType == R_PROBLEM_STRESS_MODAL || this->computeInitialAcceleration);
    for (uint m=0;m<rElement.size();m++)
    {
        for (uint i=0;i<dims;i++)
        {
            uint mp = 0;

            if (this->nodeBook.getValue(dims*rElement.getNodeId(m)+i,mp))
            {
                buffer.b[mp] += be[dims*m+i];
                if (this->computeInitialAcceleration)
                {
                    buffer.r[mp] += re[dims*m+i];
                }

                for (uint n=0;n<rElement.size();n++)
                {
                    for (uint j=0;j<dims;j++)
                    {
                        uint np = 0;

                        if (this->nodeBook.getValue(dims*rElement.getNodeId(n)+j,np))
                        {
                            buffer.A.addValue(mp,np,Ae[dims*m+i][dims*n+j]);
                            if (assembleMass)
                            {
                                // Be carries the mass in the node frames.
                                buffer.M.addValue(mp,np,Be[dims*m+i][dims*n+j]);
                            }
                            if (this->computeInitialAcceleration)
                            {
                                buffer.D.addValue(mp,np,De[dims*m+i][dims*n+j]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RSolverStress::applyLocalRotations(unsigned int elementID, RRMatrix &Ae)
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
            T.setBlock(this->localRotations[nodeId].getR(),3*i,3*i);
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

void RSolverStress::applyLocalRotations(unsigned int elementID, RRVector &fe)
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
                T.setIdentity(fe.getNRows());
                first = false;
            }
            T.setBlock(this->localRotations[nodeId].getR(),3*i,3*i);
        }
    }
    if (!first)
    {
        // T maps local to global, so with u_global = T * u_local the system
        // K_local = T^T * K_global * T has to be paired with the load vector
        // f_local = T^T * f_global. Rotating the load with T instead leaves the
        // solution satisfying equilibrium with a rotated load.
        RRMatrix Tt(T);
        Tt.transpose();

        RRVector fetmp;

        RRMatrix::mlt(Tt,fe,fetmp);
        fe = fetmp;
    }
}

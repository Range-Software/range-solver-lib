#include <atomic>
#include <cmath>

#include <omp.h>

#include "rsolveracoustic.h"
#include "rmatrixsolver.h"

RSolverAcoustic::RSolverAcoustic(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , harmonic(false)
    , angularFrequency(0.0)
{
    this->problemType = R_PROBLEM_ACOUSTICS;
}

RSolverAcoustic::~RSolverAcoustic()
{

}

bool RSolverAcoustic::hasConverged() const
{
    // The acoustic system is linear - a single solve per time step / frequency
    // is always the final answer, no outer iteration is required.
    return true;
}

double RSolverAcoustic::findNewmarkGamma()
{
    return 0.5;
}

double RSolverAcoustic::findNewmarkBeta() const
{
    // The time-march approximation coefficient is 1.0 for the backward, 0.5 for
    // the central and 0.0 for the forward approximation. Newmark requires
    // beta > 0 and beta >= 1/4 (with gamma = 1/2) for unconditional stability,
    // so the coefficient is clamped to the average-acceleration scheme.
    double beta = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient() / 2.0;
    return std::max(beta,0.25);
}

double RSolverAcoustic::findTimeStepSize() const
{
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
    if (dt < RConstants::eps)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Acoustic solver requires a positive time step size.");
    }
    return dt;
}

double RSolverAcoustic::findAverageSoundSpeed() const
{
    double cavg = 0.0;
    uint ne = 0;

    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (this->computableElements[i] && this->elementSoundSpeed[i] > RConstants::eps)
        {
            cavg += this->elementSoundSpeed[i];
            ne++;
        }
    }
    if (ne == 0)
    {
        return 0.0;
    }
    return cavg / double(ne);
}

void RSolverAcoustic::findComputableElements(RProblemType problemType)
{
    this->computableElements.resize(this->pModel->getNElements());
    this->computableElements.fill(false);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        const RMaterial &rMaterial = pElementGroup->getMaterial();

        // Density plus either an explicit speed of sound or a modulus of
        // elasticity is enough to make an element computable.
        bool hasDensity = (rMaterial.findPosition(RMaterialProperty::Density) != rMaterial.size());
        bool hasSoundSpeed = (rMaterial.findPosition(RMaterialProperty::SoundSpeed) != rMaterial.size());
        bool hasElasticity = (rMaterial.findPosition(RMaterialProperty::ModulusOfElasiticity) != rMaterial.size());

        if (hasDensity && (hasSoundSpeed || hasElasticity))
        {
            for (uint j=0;j<pElementGroup->size();j++)
            {
                this->computableElements[pElementGroup->get(j)] = true;
            }
            continue;
        }

        // Boundary entities carry no material but still take part in the solve.
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &rBoundaryCondition = pElementGroup->getBoundaryCondition(j);
            if (problemType & RBoundaryCondition::getProblemTypeMask(rBoundaryCondition.getType()))
            {
                for (uint k=0;k<pElementGroup->size();k++)
                {
                    this->computableElements[pElementGroup->get(k)] = true;
                }
                break;
            }
        }
    }
}

void RSolverAcoustic::initialize()
{
    uint nn = this->pModel->getNNodes();
    uint ne = this->pModel->getNElements();

    this->nodeVelocityPotentialImag.resize(nn,0.0);
    this->nodeAcousticPressure.resize(nn,0.0);
    this->nodeAcousticPressurePhase.resize(nn,0.0);
    this->nodeSoundPressureLevel.resize(nn,0.0);

    this->elementAcousticParticleVelocity.x.resize(ne,0.0);
    this->elementAcousticParticleVelocity.y.resize(ne,0.0);
    this->elementAcousticParticleVelocity.z.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.x.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.y.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.z.resize(ne,0.0);
    this->elementAcousticIntensity.x.resize(ne,0.0);
    this->elementAcousticIntensity.y.resize(ne,0.0);
    this->elementAcousticIntensity.z.resize(ne,0.0);
}

void RSolverAcoustic::updateScales()
{
    this->scales.setMetre(this->findMeshScale());
}

void RSolverAcoustic::recover()
{
    this->recoverVariable(R_VARIABLE_POTENTIAL,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeVelocityPotential,
                          0.0);
    this->recoverVariable(R_VARIABLE_POTENTIAL_VELOCITY,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeVelocityPotentialVelocity,
                          0.0);
    this->recoverVariable(R_VARIABLE_POTENTIAL_ACCELERATION,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeVelocityPotentialAcceleration,
                          0.0);
}

void RSolverAcoustic::generateMaterialVectors()
{
    uint ne = this->pModel->getNElements();

    RRVector elementElasticityModulus;

    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);
    this->generateMaterialVecor(RMaterialProperty::SoundSpeed,this->elementSoundSpeed);
    this->generateMaterialVecor(RMaterialProperty::ModulusOfElasiticity,elementElasticityModulus);
    this->generateMaterialVecor(RMaterialProperty::AcousticDampingFactor,this->elementDampingFactor);

    this->elementSoundSpeed.resize(ne,0.0);
    this->elementDensity.resize(ne,0.0);
    this->elementDampingFactor.resize(ne,0.0);

    // An explicitly given speed of sound wins over the one derived from the
    // modulus of elasticity and the density.
    for (uint i=0;i<ne;i++)
    {
        if (this->elementSoundSpeed[i] > RConstants::eps)
        {
            continue;
        }
        if (elementElasticityModulus[i] > RConstants::eps && this->elementDensity[i] > RConstants::eps)
        {
            this->elementSoundSpeed[i] = std::sqrt(elementElasticityModulus[i]/this->elementDensity[i]);
        }
        else
        {
            this->elementSoundSpeed[i] = 0.0;
        }
    }

    // Boundary entities carry no material of their own but the impedance and
    // radiation terms need the properties of the fluid they are bounding, so
    // fall back to the domain average.
    double cAvg = this->findAverageSoundSpeed();
    double roAvg = 0.0;
    uint nRo = 0;
    for (uint i=0;i<ne;i++)
    {
        if (this->computableElements[i] && this->elementDensity[i] > RConstants::eps)
        {
            roAvg += this->elementDensity[i];
            nRo++;
        }
    }
    if (nRo > 0)
    {
        roAvg /= double(nRo);
    }

    if (cAvg < RConstants::eps)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "No element with a valid speed of sound was found. Assign a density and either a speed of sound or a modulus of elasticity.");
    }

    for (uint i=0;i<ne;i++)
    {
        if (this->elementSoundSpeed[i] < RConstants::eps)
        {
            this->elementSoundSpeed[i] = cAvg;
        }
        if (this->elementDensity[i] < RConstants::eps)
        {
            this->elementDensity[i] = roAvg;
        }
    }
}

void RSolverAcoustic::generateBoundaryConditionVector(RBoundaryConditionType boundaryConditionType,
                                                      RVariableType variableType,
                                                      RRVector &values,
                                                      RBVector &setValues) const
{
    uint ne = this->pModel->getNElements();

    values.resize(ne,0.0);
    values.fill(0.0);
    setValues.resize(ne,false);
    setValues.fill(false);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (bc.getType() != boundaryConditionType)
            {
                continue;
            }
            uint componentPosition = bc.findComponentPosition(variableType);
            if (componentPosition == RConstants::eod)
            {
                continue;
            }
            double value = bc.getComponent(componentPosition).get(this->pModel->getTimeSolver().getCurrentTime());

            for (uint k=0;k<pElementGroup->size();k++)
            {
                values[pElementGroup->get(k)] = value;
                setValues[pElementGroup->get(k)] = true;
            }
        }
    }
}

void RSolverAcoustic::generateBoundaryDamping()
{
    uint ne = this->pModel->getNElements();

    RRVector elementAbsorption;
    RBVector elementAbsorptionSet;
    RRVector elementImpedance;
    RBVector elementImpedanceSet;

    this->generateBoundaryConditionVector(R_BOUNDARY_CONDITION_ABSORBING_BOUNDARY,
                                          R_VARIABLE_ACOUSTIC_ABSORPTION_COEFFICIENT,
                                          elementAbsorption,
                                          elementAbsorptionSet);
    this->generateBoundaryConditionVector(R_BOUNDARY_CONDITION_ACOUSTIC_IMPEDANCE,
                                          R_VARIABLE_ACOUSTIC_IMPEDANCE,
                                          elementImpedance,
                                          elementImpedanceSet);

    this->elementBoundaryDamping.resize(ne,0.0);
    this->elementBoundaryDamping.fill(0.0);

    for (uint i=0;i<ne;i++)
    {
        if (!this->computableElements[i])
        {
            continue;
        }

        double c = this->elementSoundSpeed[i];
        double ro = this->elementDensity[i];

        // Specific acoustic impedance ratio zeta = Z / (rho * c). zeta = 1 is a
        // perfectly absorbing (anechoic) boundary, zeta -> infinity is a rigid
        // wall which contributes no damping at all.
        double zeta = 0.0;

        if (elementImpedanceSet[i])
        {
            if (ro > RConstants::eps && c > RConstants::eps)
            {
                zeta = elementImpedance[i] / (ro * c);
            }
        }
        else if (elementAbsorptionSet[i])
        {
            double alpha = std::min(std::max(elementAbsorption[i],0.0),1.0);
            double r = std::sqrt(1.0 - alpha);
            if (r < 1.0 - RConstants::eps)
            {
                zeta = (1.0 + r) / (1.0 - r);
            }
        }

        if (zeta > RConstants::eps)
        {
            // C = c^2 * rho / Z = c / zeta  per unit boundary area.
            this->elementBoundaryDamping[i] = c / zeta;
        }
    }
}

void RSolverAcoustic::generateElementVectors()
{
    RBVector elementVelocityNormalSet;

    this->generateMaterialVectors();
    this->generateBoundaryDamping();
    this->generateBoundaryConditionVector(R_BOUNDARY_CONDITION_VELOCITY_NORMAL,
                                          R_VARIABLE_VELOCITY,
                                          this->elementVelocityNormal,
                                          elementVelocityNormalSet);
}

void RSolverAcoustic::prepare()
{
    const RAcousticSetup &acousticSetup = this->pModel->getProblemSetup().getAcousticSetup();

    this->harmonic = (acousticSetup.getAnalysisType() == R_ACOUSTIC_ANALYSIS_HARMONIC);
    this->angularFrequency = this->harmonic ? acousticSetup.getAngularFrequency() : 0.0;

    if (this->harmonic)
    {
        if (this->angularFrequency < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Harmonic acoustic analysis requires a positive frequency.");
        }
    }
    else if (!this->pModel->getTimeSolver().getEnabled())
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Transient acoustic analysis requires an enabled time solver. "
                     "Enable time stepping or switch the acoustic setup to a harmonic analysis.");
    }

    // Remember the state at the beginning of the step before the boundary
    // conditions of the current step overwrite it - the Newmark predictor needs
    // the previous value, also on prescribed nodes.
    this->nodeVelocityPotentialOld = this->nodeVelocityPotential;

    this->generateNodeBook(R_PROBLEM_ACOUSTICS);
    this->generateElementVectors();

    // Prescribed velocity potential (boundary and initial conditions).
    RRVector elementVelocityPotential;
    RBVector velocityPotentialSetValues;

    this->generateVariableVector(R_VARIABLE_POTENTIAL,elementVelocityPotential,velocityPotentialSetValues,true,this->firstRun,this->firstRun);
    this->pModel->convertElementToNodeVector(elementVelocityPotential,velocityPotentialSetValues,this->nodeVelocityPotential,true);

    if (this->harmonic)
    {
        // Each frequency is solved independently - the imaginary part starts
        // from zero and only prescribed nodes deviate from it.
        this->nodeVelocityPotentialImag.resize(this->pModel->getNNodes(),0.0);
        this->nodeVelocityPotentialImag.fill(0.0);
    }
    else if (this->firstRun)
    {
        // Initial condition for the first time derivative of the potential -
        // this is what puts a non-zero initial pressure into the domain.
        RRVector elementVelocityPotentialVelocity;
        RBVector velocityPotentialVelocitySetValues;

        this->generateVariableVector(R_VARIABLE_POTENTIAL_VELOCITY,
                                     elementVelocityPotentialVelocity,
                                     velocityPotentialVelocitySetValues,
                                     false,
                                     true,
                                     false);
        this->pModel->convertElementToNodeVector(elementVelocityPotentialVelocity,
                                                 velocityPotentialVelocitySetValues,
                                                 this->nodeVelocityPotentialVelocity,
                                                 true);
    }

    uint nEnabled = this->nodeBook.getNEnabled();
    uint nSystem = this->harmonic ? 2 * nEnabled : nEnabled;

    this->b.resize(nSystem);
    this->x.resize(nSystem);

    this->A.clear();
    this->A.setNRows(nSystem);
    this->b.fill(0.0);
    this->x.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b once at the end.
    int np = omp_get_max_threads();
    std::vector<RSparseMatrix> Ap(np);
    std::vector<RRVector> bp(np);
    for (int t=0;t<np;t++)
    {
        Ap[t].setNRows(nSystem);
        bp[t].resize(nSystem);
        bp[t].fill(0.0);
    }

    bool timeSolverEnabled = this->harmonic || this->pModel->getTimeSolver().getEnabled();

    // Prepare point elements.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        RPoint &point = this->pModel->getPoint(i);

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
                uint elementID = point.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_POINT(element.getType()));
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ce(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());

                Me.fill(0.0);
                Ce.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                double c = this->elementSoundSpeed[elementID];
                double beta = this->elementDampingFactor[elementID];
                double boundaryDamping = this->elementBoundaryDamping[elementID];

                // A point element has no integration points - its contribution
                // is lumped onto its single node. The assigned point volume acts
                // as the measure of the domain terms, boundary terms act over a
                // unit area, exactly as they do for line and surface elements.
                if (timeSolverEnabled)
                {
                    Me[0][0] += point.getVolume();
                }
                Ce[0][0] += beta * point.getVolume() + boundaryDamping;
                fe[0] += c * c * this->elementVelocityNormal[elementID];

                this->assemblyMatrix(elementID,Me,Ce,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
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
                RRMatrix Ce(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),1);

                Me.fill(0.0);
                Ce.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                double c = this->elementSoundSpeed[elementID];
                double beta = this->elementDampingFactor[elementID];
                double boundaryDamping = this->elementBoundaryDamping[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += dN[m][0]*J[0][0];
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<element.size();n++)
                        {
                            double NmNn = N[m] * N[n] * dV;

                            // Stiffness
                            Ke[m][n] += (B[m][0]*B[n][0]) * line.getCrossArea() * c * c * dV;

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += NmNn * line.getCrossArea();
                            }
                            // Damping - bulk absorption and boundary impedance.
                            Ce[m][n] += beta * NmNn * line.getCrossArea() + boundaryDamping * NmNn;
                        }
                        // Prescribed normal velocity source.
                        fe[m] += c * c * this->elementVelocityNormal[elementID] * N[m] * dV;
                    }
                }
                this->assemblyMatrix(elementID,Me,Ce,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
                uint elementID = surface.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ce(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),2);

                Me.fill(0.0);
                Ce.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                double c = this->elementSoundSpeed[elementID];
                double beta = this->elementDampingFactor[elementID];
                double boundaryDamping = this->elementBoundaryDamping[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<element.size();n++)
                        {
                            double NmNn = N[m] * N[n] * dV;

                            // Stiffness
                            Ke[m][n] += (B[m][0]*B[n][0]+B[m][1]*B[n][1]) * surface.getThickness() * c * c * dV;

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += NmNn * surface.getThickness();
                            }
                            // Damping - bulk absorption and boundary impedance.
                            Ce[m][n] += beta * NmNn * surface.getThickness() + boundaryDamping * NmNn;
                        }
                        // Prescribed normal velocity source.
                        fe[m] += c * c * this->elementVelocityNormal[elementID] * N[m] * dV;
                    }
                }
                this->assemblyMatrix(elementID,Me,Ce,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
                RRMatrix Ce(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),3);

                Me.fill(0.0);
                Ce.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                double c = this->elementSoundSpeed[elementID];
                double beta = this->elementDampingFactor[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<element.size();n++)
                        {
                            double NmNn = N[m] * N[n] * dV;

                            // Stiffness
                            Ke[m][n] += (B[m][0]*B[n][0] + B[m][1]*B[n][1] + B[m][2]*B[n][2]) * c * c * dV;

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += NmNn;
                            }
                            // Damping - bulk absorption.
                            Ce[m][n] += beta * NmNn;
                        }
                    }
                }
                this->assemblyMatrix(elementID,Me,Ce,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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

void RSolverAcoustic::assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    if (this->harmonic)
    {
        this->assemblyMatrixHarmonic(elementID,Me,Ce,Ke,fe,Ap,bp);
    }
    else
    {
        this->assemblyMatrixTransient(elementID,Me,Ce,Ke,fe,Ap,bp);
    }
}

void RSolverAcoustic::assemblyMatrixTransient(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    double gamma = RSolverAcoustic::findNewmarkGamma();
    double beta = this->findNewmarkBeta();
    double dt = this->findTimeStepSize();

    const RElement &element = this->pModel->getElement(elementID);

    RRMatrix Ae(element.size(),element.size());
    RRVector be(element.size());

    Ae.fill(0.0);
    be.fill(0.0);

    double a0 = 1.0 / (beta * std::pow(dt,2));
    double a1 = gamma / (beta * dt);
    double a2 = 1.0 / (beta * dt);
    double a3 = (1.0 / (2.0 * beta)) - 1.0;
    double a4 = (gamma / beta) - 1.0;
    double a5 = (1.0 / 2.0) * dt * ((gamma / beta) - 2.0);

    for (uint m=0;m<element.size();m++)
    {
        be[m] = fe[m];
        for (uint n=0;n<element.size();n++)
        {
            // The Newmark predictor uses the state at the beginning of the step.
            double pu = this->nodeVelocityPotentialOld[element.getNodeId(n)];
            double pv = this->nodeVelocityPotentialVelocity[element.getNodeId(n)];
            double pa = this->nodeVelocityPotentialAcceleration[element.getNodeId(n)];

            Ae[m][n] = Ke[m][n] + a0 * Me[m][n] + a1 * Ce[m][n];
            be[m] += Me[m][n] * (a0 * pu + a2 * pv + a3 * pa) + Ce[m][n] * (a1 * pu + a4 * pv + a5 * pa);
        }
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
                be[n] -= Ae[n][m] * this->nodeVelocityPotential[nodeID];
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

void RSolverAcoustic::assemblyMatrixHarmonic(uint elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    // (K - w^2 * M + i*w*C) * (Phi_r + i*Phi_i) = F_r + i*F_i is solved as
    //
    //     [  S   -w*C ] [Phi_r]   [F_r]
    //     [ w*C    S  ] [Phi_i] = [F_i]      with S = K - w^2 * M
    double w = this->angularFrequency;
    uint nEnabled = this->nodeBook.getNEnabled();

    const RElement &element = this->pModel->getElement(elementID);

    RRMatrix Se(element.size(),element.size());
    RRMatrix De(element.size(),element.size());
    RRVector ber(element.size());
    RRVector bei(element.size());

    Se.fill(0.0);
    De.fill(0.0);
    ber.fill(0.0);
    bei.fill(0.0);

    for (uint m=0;m<element.size();m++)
    {
        ber[m] = fe[m];
        for (uint n=0;n<element.size();n++)
        {
            Se[m][n] = Ke[m][n] - w * w * Me[m][n];
            De[m][n] = w * Ce[m][n];
        }
    }

    // Apply explicit boundary conditions.
    for (uint m=0;m<element.size();m++)
    {
        uint position;
        uint nodeID = element.getNodeId(m);
        if (!this->nodeBook.getValue(nodeID,position))
        {
            double vr = this->nodeVelocityPotential[nodeID];
            double vi = this->nodeVelocityPotentialImag[nodeID];

            for (uint n=0;n<element.size();n++)
            {
                ber[n] -= Se[n][m] * vr - De[n][m] * vi;
                bei[n] -= De[n][m] * vr + Se[n][m] * vi;
            }
        }
    }

    // Assembly final matrix system
    for (uint m=0;m<element.size();m++)
    {
        uint mp;

        if (this->nodeBook.getValue(element.getNodeId(m),mp))
        {
            bp[mp] += ber[m];
            bp[mp+nEnabled] += bei[m];

            for (uint n=0;n<element.size();n++)
            {
                uint np = 0;

                if (this->nodeBook.getValue(element.getNodeId(n),np))
                {
                    Ap.addValue(mp,np,Se[m][n]);
                    Ap.addValue(mp,np+nEnabled,-De[m][n]);
                    Ap.addValue(mp+nEnabled,np,De[m][n]);
                    Ap.addValue(mp+nEnabled,np+nEnabled,Se[m][n]);
                }
            }
        }
    }
}

void RSolverAcoustic::solve()
{
    if (this->harmonic)
    {
        this->solveHarmonic();
    }
    else
    {
        this->solveTransient();
    }
}

void RSolverAcoustic::solveTransient()
{
    try
    {
        RLogger::indent();
        // K + a0*M + a1*C is symmetric positive definite.
        RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,1);
        RLogger::unindent();
    }
    catch (const RError &)
    {
        RLogger::unindent();
        throw;
    }

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        if (this->nodeBook.getValue(i,position))
        {
            this->nodeVelocityPotential[i] = this->x[position];
        }
    }

    double gamma = RSolverAcoustic::findNewmarkGamma();
    double beta = this->findNewmarkBeta();
    double dt = this->findTimeStepSize();

    double a0 = 1.0 / (beta * std::pow(dt,2));
    double a2 = 1.0 / (beta * dt);
    double a3 = (1.0 / (2.0 * beta)) - 1.0;
    double a6 = dt * (1.0 - gamma);
    double a7 = dt * gamma;

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        double puo = this->nodeVelocityPotentialOld[i];
        double pvo = this->nodeVelocityPotentialVelocity[i];
        double pao = this->nodeVelocityPotentialAcceleration[i];

        double pu = this->nodeVelocityPotential[i];
        double pa = a0 * (pu - puo) - a2 * pvo - a3 * pao;
        double pv = pvo + a6 * pao + a7 * pa;

        this->nodeVelocityPotentialVelocity[i] = pv;
        this->nodeVelocityPotentialAcceleration[i] = pa;
    }
}

void RSolverAcoustic::solveHarmonic()
{
    // The real block form of the complex system is not symmetric and, away from
    // the static limit, indefinite. Restarted GMRES stagnates on it whenever the
    // Krylov space is much smaller than the system, so widen the restart as far
    // as a fixed memory budget for the two Krylov bases allows. The configured
    // value is only ever increased, never reduced.
    RMatrixSolverConf matrixSolverConf = this->pModel->getMatrixSolverConf(RMatrixSolverConf::GMRES);

    uint nSystem = uint(this->b.size());
    if (nSystem > 0)
    {
        const double krylovMemoryBudget = 256.0 * 1024.0 * 1024.0;
        uint maxInnerIterations = uint(std::max(1.0,krylovMemoryBudget / (2.0 * double(sizeof(double)) * double(nSystem))));
        uint nInnerIterations = std::min(nSystem+1,maxInnerIterations);

        if (matrixSolverConf.getNInnerIterations() < nInnerIterations)
        {
            RLogger::info("Increasing GMRES restart from %u to %u for the frequency domain system.\n",
                          matrixSolverConf.getNInnerIterations(),
                          nInnerIterations);
            matrixSolverConf.setNInnerIterations(nInnerIterations);
        }
    }

    try
    {
        RLogger::indent();
        RMatrixSolver matrixSolver(matrixSolverConf);
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,1);
        RLogger::unindent();
    }
    catch (const RError &)
    {
        RLogger::unindent();
        throw;
    }

    uint nEnabled = this->nodeBook.getNEnabled();

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        if (this->nodeBook.getValue(i,position))
        {
            this->nodeVelocityPotential[i] = this->x[position];
            this->nodeVelocityPotentialImag[i] = this->x[position+nEnabled];
        }
    }

    // Keep the time derivatives consistent with the harmonic solution so that
    // monitoring points and derived quantities stay meaningful.
    double w = this->angularFrequency;
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        double pr = this->nodeVelocityPotential[i];
        double pi = this->nodeVelocityPotentialImag[i];

        // d/dt -> i*w  and  d2/dt2 -> -w^2 (amplitudes).
        this->nodeVelocityPotentialVelocity[i] = w * std::sqrt(pr*pr + pi*pi);
        this->nodeVelocityPotentialAcceleration[i] = w * w * std::sqrt(pr*pr + pi*pi);
    }
}

void RSolverAcoustic::process()
{
    this->processAcousticPressure();
    this->processAcousticParticleVelocity();
    this->processAcousticIntensity();
}

void RSolverAcoustic::processAcousticPressure()
{
    uint nn = this->pModel->getNNodes();

    this->nodeAcousticPressure.resize(nn,0.0);
    this->nodeAcousticPressurePhase.resize(nn,0.0);
    this->nodeSoundPressureLevel.resize(nn,0.0);

    // No element value is prescribed - let the conversion produce a distance
    // weighted average of the surrounding elements.
    RRVector nodeDensity;
    RBVector nodeDensitySetValues(this->pModel->getNElements(),false);
    this->pModel->convertElementToNodeVector(this->elementDensity,nodeDensitySetValues,nodeDensity);
    nodeDensity.resize(nn,0.0);

    // Nodes which are not touched by any element carrying a density (e.g. nodes
    // of a zero thickness boundary surface) fall back to the domain average.
    double roAvg = 0.0;
    uint nRo = 0;
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (this->computableElements[i] && this->elementDensity[i] > RConstants::eps)
        {
            roAvg += this->elementDensity[i];
            nRo++;
        }
    }
    if (nRo > 0)
    {
        roAvg /= double(nRo);
        for (uint i=0;i<nn;i++)
        {
            if (nodeDensity[i] < RConstants::eps)
            {
                nodeDensity[i] = roAvg;
            }
        }
    }

    double pRef = this->pModel->getProblemSetup().getAcousticSetup().getReferencePressure();
    if (pRef < RConstants::eps)
    {
        pRef = R_ACOUSTIC_REFERENCE_PRESSURE;
    }

    double w = this->angularFrequency;

    for (uint i=0;i<nn;i++)
    {
        double ro = nodeDensity[i];
        double pAmplitude = 0.0;

        if (this->harmonic)
        {
            // p = rho * dphi/dt -> P = i*w*rho*Phi
            double pr = -w * ro * this->nodeVelocityPotentialImag[i];
            double pi =  w * ro * this->nodeVelocityPotential[i];

            this->nodeAcousticPressure[i] = std::sqrt(pr*pr + pi*pi);
            this->nodeAcousticPressurePhase[i] = std::atan2(pi,pr) * 180.0 / RConstants::pi;
            // The sound pressure level of a harmonic signal is based on the RMS
            // value, which is the amplitude divided by sqrt(2).
            pAmplitude = this->nodeAcousticPressure[i] / std::sqrt(2.0);
        }
        else
        {
            // p = rho * dphi/dt
            this->nodeAcousticPressure[i] = ro * this->nodeVelocityPotentialVelocity[i];
            this->nodeAcousticPressurePhase[i] = 0.0;
            pAmplitude = std::fabs(this->nodeAcousticPressure[i]);
        }

        this->nodeSoundPressureLevel[i] = (pAmplitude > RConstants::eps)
                                        ? 20.0 * std::log10(pAmplitude/pRef)
                                        : 0.0;
    }
}

void RSolverAcoustic::processAcousticParticleVelocity()
{
    uint ne = this->pModel->getNElements();

    this->elementAcousticParticleVelocity.x.resize(ne,0.0);
    this->elementAcousticParticleVelocity.y.resize(ne,0.0);
    this->elementAcousticParticleVelocity.z.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.x.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.y.resize(ne,0.0);
    this->elementAcousticParticleVelocityImag.z.resize(ne,0.0);

    this->elementAcousticParticleVelocity.x.fill(0.0);
    this->elementAcousticParticleVelocity.y.fill(0.0);
    this->elementAcousticParticleVelocity.z.fill(0.0);
    this->elementAcousticParticleVelocityImag.x.fill(0.0);
    this->elementAcousticParticleVelocityImag.y.fill(0.0);
    this->elementAcousticParticleVelocityImag.z.fill(0.0);

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

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
                uint elementID = line.get(j);

                if (!this->computableElements[elementID] || line.getCrossArea() < RConstants::eps)
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRVector B(element.size());

                B.fill(0.0);
                double volume = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m] += dN[m][0] * J[0][0] * dV;
                    }
                    volume += dV;
                }

                if (volume < RConstants::eps)
                {
                    continue;
                }
                for (uint m=0;m<B.size();m++)
                {
                    B[m] /= volume;
                }

                double vir = 0.0;
                double vii = 0.0;

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    // u = -grad(phi)
                    vir -= B[k] * this->nodeVelocityPotential[nodeID];
                    vii -= B[k] * this->nodeVelocityPotentialImag[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                this->elementAcousticParticleVelocity.x[elementID] = R[0][0]*vir;
                this->elementAcousticParticleVelocity.y[elementID] = R[1][0]*vir;
                this->elementAcousticParticleVelocity.z[elementID] = R[2][0]*vir;
                this->elementAcousticParticleVelocityImag.x[elementID] = R[0][0]*vii;
                this->elementAcousticParticleVelocityImag.y[elementID] = R[1][0]*vii;
                this->elementAcousticParticleVelocityImag.z[elementID] = R[2][0]*vii;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process acoustic particle velocity.");
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

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
                uint elementID = surface.get(j);

                if (!this->computableElements[elementID] || surface.getThickness() < RConstants::eps)
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix B(element.size(),2);

                B.fill(0.0);
                double volume = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]) * dV;
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]) * dV;
                    }
                    volume += dV;
                }

                if (volume < RConstants::eps)
                {
                    continue;
                }
                for (uint m=0;m<B.getNRows();m++)
                {
                    B[m][0] /= volume;
                    B[m][1] /= volume;
                }

                double vir = 0.0;
                double vjr = 0.0;
                double vii = 0.0;
                double vji = 0.0;

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    // u = -grad(phi)
                    vir -= B[k][0] * this->nodeVelocityPotential[nodeID];
                    vjr -= B[k][1] * this->nodeVelocityPotential[nodeID];
                    vii -= B[k][0] * this->nodeVelocityPotentialImag[nodeID];
                    vji -= B[k][1] * this->nodeVelocityPotentialImag[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                this->elementAcousticParticleVelocity.x[elementID] = R[0][0]*vir + R[0][1]*vjr;
                this->elementAcousticParticleVelocity.y[elementID] = R[1][0]*vir + R[1][1]*vjr;
                this->elementAcousticParticleVelocity.z[elementID] = R[2][0]*vir + R[2][1]*vjr;
                this->elementAcousticParticleVelocityImag.x[elementID] = R[0][0]*vii + R[0][1]*vji;
                this->elementAcousticParticleVelocityImag.y[elementID] = R[1][0]*vii + R[1][1]*vji;
                this->elementAcousticParticleVelocityImag.z[elementID] = R[2][0]*vii + R[2][1]*vji;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process acoustic particle velocity.");
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
                uint elementID = volume.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix B(element.size(),3);

                B.fill(0.0);
                double elementVolume = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]) * dV;
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]) * dV;
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]) * dV;
                    }
                    elementVolume += dV;
                }

                if (elementVolume < RConstants::eps)
                {
                    continue;
                }
                for (uint m=0;m<B.getNRows();m++)
                {
                    B[m][0] /= elementVolume;
                    B[m][1] /= elementVolume;
                    B[m][2] /= elementVolume;
                }

                RR3Vector ver(0.0,0.0,0.0);
                RR3Vector vei(0.0,0.0,0.0);

                for (uint m=0;m<element.size();m++)
                {
                    uint nodeId = element.getNodeId(m);

                    // u = -grad(phi)
                    ver[0] -= B[m][0] * this->nodeVelocityPotential[nodeId];
                    ver[1] -= B[m][1] * this->nodeVelocityPotential[nodeId];
                    ver[2] -= B[m][2] * this->nodeVelocityPotential[nodeId];
                    vei[0] -= B[m][0] * this->nodeVelocityPotentialImag[nodeId];
                    vei[1] -= B[m][1] * this->nodeVelocityPotentialImag[nodeId];
                    vei[2] -= B[m][2] * this->nodeVelocityPotentialImag[nodeId];
                }

                this->elementAcousticParticleVelocity.x[elementID] = ver[0];
                this->elementAcousticParticleVelocity.y[elementID] = ver[1];
                this->elementAcousticParticleVelocity.z[elementID] = ver[2];
                this->elementAcousticParticleVelocityImag.x[elementID] = vei[0];
                this->elementAcousticParticleVelocityImag.y[elementID] = vei[1];
                this->elementAcousticParticleVelocityImag.z[elementID] = vei[2];
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process acoustic particle velocity.");
        }
    }
}

void RSolverAcoustic::processAcousticIntensity()
{
    uint ne = this->pModel->getNElements();

    this->elementAcousticIntensity.x.resize(ne,0.0);
    this->elementAcousticIntensity.y.resize(ne,0.0);
    this->elementAcousticIntensity.z.resize(ne,0.0);

    this->elementAcousticIntensity.x.fill(0.0);
    this->elementAcousticIntensity.y.fill(0.0);
    this->elementAcousticIntensity.z.fill(0.0);

    double w = this->angularFrequency;

    for (uint i=0;i<ne;i++)
    {
        if (!this->computableElements[i])
        {
            continue;
        }

        const RElement &element = this->pModel->getElement(i);
        if (element.size() == 0)
        {
            continue;
        }

        // Element averaged potential.
        double phir = 0.0;
        double phii = 0.0;
        for (uint j=0;j<element.size();j++)
        {
            phir += this->nodeVelocityPotential[element.getNodeId(j)];
            phii += this->nodeVelocityPotentialImag[element.getNodeId(j)];
        }
        phir /= double(element.size());
        phii /= double(element.size());

        double ro = this->elementDensity[i];

        if (this->harmonic)
        {
            // Time averaged intensity I = 0.5 * Re{ P * conj(U) }.
            double pr = -w * ro * phii;
            double pi =  w * ro * phir;

            this->elementAcousticIntensity.x[i] = 0.5 * (pr * this->elementAcousticParticleVelocity.x[i]
                                                       + pi * this->elementAcousticParticleVelocityImag.x[i]);
            this->elementAcousticIntensity.y[i] = 0.5 * (pr * this->elementAcousticParticleVelocity.y[i]
                                                       + pi * this->elementAcousticParticleVelocityImag.y[i]);
            this->elementAcousticIntensity.z[i] = 0.5 * (pr * this->elementAcousticParticleVelocity.z[i]
                                                       + pi * this->elementAcousticParticleVelocityImag.z[i]);
        }
        else
        {
            // Instantaneous intensity I = p * u.
            double pv = 0.0;
            for (uint j=0;j<element.size();j++)
            {
                pv += this->nodeVelocityPotentialVelocity[element.getNodeId(j)];
            }
            pv /= double(element.size());

            double p = ro * pv;

            this->elementAcousticIntensity.x[i] = p * this->elementAcousticParticleVelocity.x[i];
            this->elementAcousticIntensity.y[i] = p * this->elementAcousticParticleVelocity.y[i];
            this->elementAcousticIntensity.z[i] = p * this->elementAcousticParticleVelocity.z[i];
        }
    }

    if (this->harmonic)
    {
        // The stored particle velocity is the amplitude of the complex value.
        for (uint i=0;i<ne;i++)
        {
            this->elementAcousticParticleVelocity.x[i] = std::sqrt(std::pow(this->elementAcousticParticleVelocity.x[i],2)
                                                                 + std::pow(this->elementAcousticParticleVelocityImag.x[i],2));
            this->elementAcousticParticleVelocity.y[i] = std::sqrt(std::pow(this->elementAcousticParticleVelocity.y[i],2)
                                                                 + std::pow(this->elementAcousticParticleVelocityImag.y[i],2));
            this->elementAcousticParticleVelocity.z[i] = std::sqrt(std::pow(this->elementAcousticParticleVelocity.z[i],2)
                                                                 + std::pow(this->elementAcousticParticleVelocityImag.z[i],2));
        }
    }
}

void RSolverAcoustic::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    uint nn = this->pModel->getNNodes();
    uint ne = this->pModel->getNElements();

    // Velocity potential
    uint velocityPotentialPos = this->pModel->findVariable(R_VARIABLE_POTENTIAL);
    if (velocityPotentialPos == RConstants::eod)
    {
        velocityPotentialPos = this->pModel->addVariable(R_VARIABLE_POTENTIAL);
        this->pModel->getVariable(velocityPotentialPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeVelocityPotential),
                    RStatistics::findMaximumValue(this->nodeVelocityPotential));
    }
    RVariable &velocityPotential = this->pModel->getVariable(velocityPotentialPos);

    velocityPotential.setApplyType(R_VARIABLE_APPLY_NODE);
    velocityPotential.resize(1,nn);
    for (uint i=0;i<nn;i++)
    {
        velocityPotential.setValue(0,i,this->nodeVelocityPotential[i]);
    }

    // Velocity potential - first time derivative. Also carries the Newmark
    // state over to the next time step.
    uint velocityPotentialVelocityPos = this->pModel->findVariable(R_VARIABLE_POTENTIAL_VELOCITY);
    if (velocityPotentialVelocityPos == RConstants::eod)
    {
        velocityPotentialVelocityPos = this->pModel->addVariable(R_VARIABLE_POTENTIAL_VELOCITY);
        this->pModel->getVariable(velocityPotentialVelocityPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeVelocityPotentialVelocity),
                    RStatistics::findMaximumValue(this->nodeVelocityPotentialVelocity));
    }
    RVariable &velocityPotentialVelocity = this->pModel->getVariable(velocityPotentialVelocityPos);

    velocityPotentialVelocity.setApplyType(R_VARIABLE_APPLY_NODE);
    velocityPotentialVelocity.resize(1,nn);
    for (uint i=0;i<nn;i++)
    {
        velocityPotentialVelocity.setValue(0,i,this->nodeVelocityPotentialVelocity[i]);
    }

    // Velocity potential - second time derivative.
    uint velocityPotentialAccelerationPos = this->pModel->findVariable(R_VARIABLE_POTENTIAL_ACCELERATION);
    if (velocityPotentialAccelerationPos == RConstants::eod)
    {
        velocityPotentialAccelerationPos = this->pModel->addVariable(R_VARIABLE_POTENTIAL_ACCELERATION);
        this->pModel->getVariable(velocityPotentialAccelerationPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeVelocityPotentialAcceleration),
                    RStatistics::findMaximumValue(this->nodeVelocityPotentialAcceleration));
    }
    RVariable &velocityPotentialAcceleration = this->pModel->getVariable(velocityPotentialAccelerationPos);

    velocityPotentialAcceleration.setApplyType(R_VARIABLE_APPLY_NODE);
    velocityPotentialAcceleration.resize(1,nn);
    for (uint i=0;i<nn;i++)
    {
        velocityPotentialAcceleration.setValue(0,i,this->nodeVelocityPotentialAcceleration[i]);
    }

    if (this->harmonic)
    {
        // Velocity potential - imaginary part.
        uint velocityPotentialImagPos = this->pModel->findVariable(R_VARIABLE_POTENTIAL_IMAGINARY);
        if (velocityPotentialImagPos == RConstants::eod)
        {
            velocityPotentialImagPos = this->pModel->addVariable(R_VARIABLE_POTENTIAL_IMAGINARY);
            this->pModel->getVariable(velocityPotentialImagPos).getVariableData().setMinMaxDisplayValue(
                        RStatistics::findMinimumValue(this->nodeVelocityPotentialImag),
                        RStatistics::findMaximumValue(this->nodeVelocityPotentialImag));
        }
        RVariable &velocityPotentialImag = this->pModel->getVariable(velocityPotentialImagPos);

        velocityPotentialImag.setApplyType(R_VARIABLE_APPLY_NODE);
        velocityPotentialImag.resize(1,nn);
        for (uint i=0;i<nn;i++)
        {
            velocityPotentialImag.setValue(0,i,this->nodeVelocityPotentialImag[i]);
        }

        // Acoustic phase.
        uint acousticPhasePos = this->pModel->findVariable(R_VARIABLE_ACOUSTIC_PHASE);
        if (acousticPhasePos == RConstants::eod)
        {
            acousticPhasePos = this->pModel->addVariable(R_VARIABLE_ACOUSTIC_PHASE);
            this->pModel->getVariable(acousticPhasePos).getVariableData().setMinMaxDisplayValue(-180.0,180.0);
        }
        RVariable &acousticPhase = this->pModel->getVariable(acousticPhasePos);

        acousticPhase.setApplyType(R_VARIABLE_APPLY_NODE);
        acousticPhase.resize(1,nn);
        for (uint i=0;i<nn;i++)
        {
            acousticPhase.setValue(0,i,this->nodeAcousticPressurePhase[i]);
        }
    }

    // Acoustic pressure
    uint acousticPressurePos = this->pModel->findVariable(R_VARIABLE_ACOUSTIC_PRESSURE);
    if (acousticPressurePos == RConstants::eod)
    {
        acousticPressurePos = this->pModel->addVariable(R_VARIABLE_ACOUSTIC_PRESSURE);
        this->pModel->getVariable(acousticPressurePos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeAcousticPressure),
                    RStatistics::findMaximumValue(this->nodeAcousticPressure));
    }
    RVariable &acousticPressure = this->pModel->getVariable(acousticPressurePos);

    acousticPressure.setApplyType(R_VARIABLE_APPLY_NODE);
    acousticPressure.resize(1,nn);
    for (uint i=0;i<nn;i++)
    {
        acousticPressure.setValue(0,i,this->nodeAcousticPressure[i]);
    }

    // Sound pressure level
    uint soundPressureLevelPos = this->pModel->findVariable(R_VARIABLE_ACOUSTIC_SOUND_PRESSURE_LEVEL);
    if (soundPressureLevelPos == RConstants::eod)
    {
        soundPressureLevelPos = this->pModel->addVariable(R_VARIABLE_ACOUSTIC_SOUND_PRESSURE_LEVEL);
        this->pModel->getVariable(soundPressureLevelPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeSoundPressureLevel),
                    RStatistics::findMaximumValue(this->nodeSoundPressureLevel));
    }
    RVariable &soundPressureLevel = this->pModel->getVariable(soundPressureLevelPos);

    soundPressureLevel.setApplyType(R_VARIABLE_APPLY_NODE);
    soundPressureLevel.resize(1,nn);
    for (uint i=0;i<nn;i++)
    {
        soundPressureLevel.setValue(0,i,this->nodeSoundPressureLevel[i]);
    }

    // Acoustic particle velocity
    uint acousticParticleVelocityPos = this->pModel->findVariable(R_VARIABLE_ACOUSTIC_PARTICLE_VELOCITY);
    if (acousticParticleVelocityPos == RConstants::eod)
    {
        acousticParticleVelocityPos = this->pModel->addVariable(R_VARIABLE_ACOUSTIC_PARTICLE_VELOCITY);

        double umin = 0.0;
        double umax = 0.0;
        for (uint i=0;i<this->elementAcousticParticleVelocity.x.size();i++)
        {
            double u = RR3Vector(this->elementAcousticParticleVelocity.x[i],
                                 this->elementAcousticParticleVelocity.y[i],
                                 this->elementAcousticParticleVelocity.z[i]).length();
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

        this->pModel->getVariable(acousticParticleVelocityPos).getVariableData().setMinMaxDisplayValue(umin,umax);
    }
    RVariable &acousticParticleVelocity = this->pModel->getVariable(acousticParticleVelocityPos);

    acousticParticleVelocity.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    acousticParticleVelocity.resize(3,ne);
    for (uint i=0;i<ne;i++)
    {
        acousticParticleVelocity.setValue(0,i,this->elementAcousticParticleVelocity.x[i]);
        acousticParticleVelocity.setValue(1,i,this->elementAcousticParticleVelocity.y[i]);
        acousticParticleVelocity.setValue(2,i,this->elementAcousticParticleVelocity.z[i]);
    }

    // Acoustic intensity
    uint acousticIntensityPos = this->pModel->findVariable(R_VARIABLE_ACOUSTIC_INTENSITY);
    if (acousticIntensityPos == RConstants::eod)
    {
        acousticIntensityPos = this->pModel->addVariable(R_VARIABLE_ACOUSTIC_INTENSITY);

        double imin = 0.0;
        double imax = 0.0;
        for (uint i=0;i<this->elementAcousticIntensity.x.size();i++)
        {
            double u = RR3Vector(this->elementAcousticIntensity.x[i],
                                 this->elementAcousticIntensity.y[i],
                                 this->elementAcousticIntensity.z[i]).length();
            if (i == 0)
            {
                imin = imax = u;
            }
            else
            {
                imin = std::min(imin,u);
                imax = std::max(imax,u);
            }
        }

        this->pModel->getVariable(acousticIntensityPos).getVariableData().setMinMaxDisplayValue(imin,imax);
    }
    RVariable &acousticIntensity = this->pModel->getVariable(acousticIntensityPos);

    acousticIntensity.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    acousticIntensity.resize(3,ne);
    for (uint i=0;i<ne;i++)
    {
        acousticIntensity.setValue(0,i,this->elementAcousticIntensity.x[i]);
        acousticIntensity.setValue(1,i,this->elementAcousticIntensity.y[i]);
        acousticIntensity.setValue(2,i,this->elementAcousticIntensity.z[i]);
    }

    RLogger::unindent();
}

void RSolverAcoustic::statistics()
{
    this->printStats(R_VARIABLE_POTENTIAL);
    this->printStats(R_VARIABLE_ACOUSTIC_PRESSURE);
    this->printStats(R_VARIABLE_ACOUSTIC_SOUND_PRESSURE_LEVEL);
    this->printStats(R_VARIABLE_ACOUSTIC_PARTICLE_VELOCITY);
    this->printStats(R_VARIABLE_ACOUSTIC_INTENSITY);
    if (this->harmonic)
    {
        this->printStats(R_VARIABLE_ACOUSTIC_PHASE);
    }
    this->processMonitoringPoints();
}

#include <atomic>

#include <omp.h>

#include "rsolverpotential.h"
#include "rmatrixsolver.h"

RSolverPotential::RSolverPotential(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , conductivityScale(1.0)
{
    this->problemType = R_PROBLEM_POTENTIAL;
}

RSolverPotential::~RSolverPotential()
{

}

bool RSolverPotential::hasConverged() const
{
    return true;
}

void RSolverPotential::initialize()
{
}

void RSolverPotential::updateScales()
{

}

void RSolverPotential::recover()
{
    this->recoverVariable(R_VARIABLE_PRESSURE,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodePressure,
                          0.0);
}

void RSolverPotential::prepare()
{
    RRVector elementPressure;
    RBVector pressureSetValues;

    this->generatePressureNodeBook();
    this->generateBoundaryVectors(elementPressure,pressureSetValues);

    if (this->nodeBook.getNEnabled() == this->nodeBook.size())
    {
        // Pure flux conditions determine the pressure up to a constant only.
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Flow through porous media requires at least one \'%s\' boundary condition.",
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT).toUtf8().constData());
    }

    // Gravity is an environment condition - boundary and initial conditions
    // carry no gravitational acceleration.
    RBVector gravitySetValues;
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_X,this->elementGravityX,gravitySetValues,false,false,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Y,this->elementGravityY,gravitySetValues,false,false,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Z,this->elementGravityZ,gravitySetValues,false,false,true);

    this->generateMaterialVecor(RMaterialProperty::PermeabilityToFluids,this->elementPermeability);
    this->generateMaterialVecor(RMaterialProperty::DynamicViscosity,this->elementViscosity);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);

    this->elementConductivity.resize(this->pModel->getNElements());
    this->conductivityScale = 0.0;
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        this->elementConductivity[i] = (this->elementViscosity[i] > 0.0) ? this->elementPermeability[i] / this->elementViscosity[i] : 0.0;
        if (this->computableElements[i])
        {
            this->conductivityScale = std::max(this->conductivityScale,this->elementConductivity[i]);
        }
    }
    // Typical conductivities (1e-15 - 1e-3) would leave the matrix entries
    // close to the matrix solver's round-off limits. Dividing every equation by
    // the same conductivity does not change the solution.
    if (this->conductivityScale <= 0.0)
    {
        this->conductivityScale = 1.0;
    }

    this->pModel->convertElementToNodeVector(elementPressure,pressureSetValues,this->nodePressure,true);

    this->b.resize(this->nodeBook.getNEnabled());
    this->x.resize(this->nodeBook.getNEnabled());

    this->A.clear();
    this->A.setNRows(this->b.size());
    this->b.fill(0.0);
    this->x.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b once at the end.
    int np = omp_get_max_threads();
    std::vector<RSparseMatrix> Ap(np);
    std::vector<RRVector> bp(np);
    for (int t=0;t<np;t++)
    {
        Ap[t].setNRows(this->b.size());
        bp[t].resize(this->b.size());
        bp[t].fill(0.0);
    }

    // Weak form of div(K*(grad(p) - rho*g)) = 0 with K = k/mu:
    //   int(grad(N).K.grad(p)) = int(grad(N).K*rho*g) + int_G(N*v_in)
    // where v_in is the normal inflow velocity on the boundary G.

    // Prepare line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        if (line.getCrossArea() <= 0.0)
        {
            continue;
        }

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
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),1);

                Ke.fill(0.0);
                fe.fill(0.0);

                double K = this->elementConductivity[elementID] / this->conductivityScale;

                // Gravity projected to the element direction.
                RRMatrix R;
                RRVector t;
                element.findTransformationMatrix(this->pModel->getNodes(),R,t);
                RR3Vector g(this->getElementGravity(elementID));
                double rhoG = this->elementDensity[elementID] * (R[0][0]*g[0] + R[1][0]*g[1] + R[2][0]*g[2]);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW() * line.getCrossArea();

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += dN[m][0]*J[0][0];
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            Ke[m][n] += B[m][0] * B[n][0] * K * dV;
                        }
                        fe[m] += B[m][0] * K * rhoG * dV;
                    }
                }
                this->assemblyMatrix(elementID,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
            try
            {
                uint elementID = surface.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                bool hasConduction = (surface.getThickness() > 0.0 && this->elementConductivity[elementID] > 0.0);
                bool hasInflow = (this->elementInflowVelocity[elementID] != 0.0);

                if (!hasConduction && !hasInflow)
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),2);

                Ke.fill(0.0);
                fe.fill(0.0);

                double K = this->elementConductivity[elementID] / this->conductivityScale;

                // Gravity projected to the element plane.
                RRMatrix R;
                RRVector t;
                element.findTransformationMatrix(this->pModel->getNodes(),R,t);
                RR3Vector g(this->getElementGravity(elementID));
                double rhoGi = this->elementDensity[elementID] * (R[0][0]*g[0] + R[1][0]*g[1] + R[2][0]*g[2]);
                double rhoGj = this->elementDensity[elementID] * (R[0][1]*g[0] + R[1][1]*g[1] + R[2][1]*g[2]);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dA = detJ * shapeFunc.getW();

                    if (hasConduction)
                    {
                        double dV = dA * surface.getThickness();

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
                                Ke[m][n] += (B[m][0] * B[n][0] + B[m][1] * B[n][1]) * K * dV;
                            }
                            fe[m] += (B[m][0] * rhoGi + B[m][1] * rhoGj) * K * dV;
                        }
                    }

                    if (hasInflow)
                    {
                        for (unsigned m=0;m<element.size();m++)
                        {
                            fe[m] += this->elementInflowVelocity[elementID] / this->conductivityScale * N[m] * dA;
                        }
                    }
                }
                this->assemblyMatrix(elementID,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),3);

                Ke.fill(0.0);
                fe.fill(0.0);

                double K = this->elementConductivity[elementID] / this->conductivityScale;
                RR3Vector rhoG(this->getElementGravity(elementID));
                rhoG *= this->elementDensity[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    double dV = detJ * shapeFunc.getW();

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
                            Ke[m][n] += (B[m][0]*B[n][0] + B[m][1]*B[n][1] + B[m][2]*B[n][2]) * K * dV;
                        }
                        fe[m] += (B[m][0]*rhoG[0] + B[m][1]*rhoG[1] + B[m][2]*rhoG[2]) * K * dV;
                    }
                }
                this->assemblyMatrix(elementID,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
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

void RSolverPotential::solve()
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

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        if (this->nodeBook.getValue(i,position))
        {
            this->nodePressure[i] = this->x[position];
        }
    }
}

void RSolverPotential::process()
{
    this->elementVelocity.assign(this->pModel->getNElements(),RR3Vector(0.0,0.0,0.0));

    // Velocity is a density - the shape function derivatives are averaged over
    // the integration points and must not be weighted by the Jacobian
    // determinant.

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        if (line.getCrossArea() <= 0.0)
        {
            continue;
        }

        for (uint j=0;j<line.size();j++)
        {
            uint elementID = line.get(j);

            if (!this->computableElements[elementID] || this->elementConductivity[elementID] <= 0.0)
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRVector B(element.size(),0.0);

            for (uint k=0;k<nInp;k++)
            {
                const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                const RRMatrix &dN = shapeFunc.getDN();
                RRMatrix J, Rt;
                element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                for (uint m=0;m<dN.getNRows();m++)
                {
                    B[m] += dN[m][0] * J[0][0] / double(nInp);
                }
            }

            RRMatrix R;
            RRVector t;
            element.findTransformationMatrix(this->pModel->getNodes(),R,t);
            RR3Vector g(this->getElementGravity(elementID));

            double vi = this->elementDensity[elementID] * (R[0][0]*g[0] + R[1][0]*g[1] + R[2][0]*g[2]);
            for (uint k=0;k<element.size();k++)
            {
                vi -= B[k] * this->nodePressure[element.getNodeId(k)];
            }
            vi *= this->elementConductivity[elementID];

            this->elementVelocity[elementID] = RR3Vector(R[0][0]*vi,R[1][0]*vi,R[2][0]*vi);
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

        if (surface.getThickness() <= 0.0)
        {
            continue;
        }

        for (uint j=0;j<surface.size();j++)
        {
            uint elementID = surface.get(j);

            if (!this->computableElements[elementID] || this->elementConductivity[elementID] <= 0.0)
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRMatrix B(element.size(),2);

            B.fill(0.0);
            for (uint k=0;k<nInp;k++)
            {
                const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                const RRMatrix &dN = shapeFunc.getDN();
                RRMatrix J, Rt;
                element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                for (uint m=0;m<dN.getNRows();m++)
                {
                    B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]) / double(nInp);
                    B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]) / double(nInp);
                }
            }

            RRMatrix R;
            RRVector t;
            element.findTransformationMatrix(this->pModel->getNodes(),R,t);
            RR3Vector g(this->getElementGravity(elementID));

            double vi = this->elementDensity[elementID] * (R[0][0]*g[0] + R[1][0]*g[1] + R[2][0]*g[2]);
            double vj = this->elementDensity[elementID] * (R[0][1]*g[0] + R[1][1]*g[1] + R[2][1]*g[2]);
            for (uint k=0;k<element.size();k++)
            {
                uint nodeID = element.getNodeId(k);
                vi -= B[k][0] * this->nodePressure[nodeID];
                vj -= B[k][1] * this->nodePressure[nodeID];
            }
            vi *= this->elementConductivity[elementID];
            vj *= this->elementConductivity[elementID];

            this->elementVelocity[elementID] = RR3Vector(R[0][0]*vi + R[0][1]*vj,
                                                         R[1][0]*vi + R[1][1]*vj,
                                                         R[2][0]*vi + R[2][1]*vj);
        }
    }

    // Process volume elements.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        RVolume &volume = this->pModel->getVolume(i);

        for (uint j=0;j<volume.size();j++)
        {
            uint elementID = volume.get(j);

            if (!this->computableElements[elementID] || this->elementConductivity[elementID] <= 0.0)
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRMatrix B(element.size(),3);

            B.fill(0.0);
            for (uint k=0;k<nInp;k++)
            {
                const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                const RRMatrix &dN = shapeFunc.getDN();
                RRMatrix J, Rt;
                element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                for (uint m=0;m<dN.getNRows();m++)
                {
                    B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]) / double(nInp);
                    B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]) / double(nInp);
                    B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]) / double(nInp);
                }
            }

            RR3Vector velocity(this->getElementGravity(elementID));
            velocity *= this->elementDensity[elementID];
            for (uint m=0;m<element.size();m++)
            {
                uint nodeID = element.getNodeId(m);
                velocity[0] -= B[m][0] * this->nodePressure[nodeID];
                velocity[1] -= B[m][1] * this->nodePressure[nodeID];
                velocity[2] -= B[m][2] * this->nodePressure[nodeID];
            }
            velocity *= this->elementConductivity[elementID];

            this->elementVelocity[elementID] = velocity;
        }
    }
}

void RSolverPotential::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Pressure
    uint pressurePos = this->pModel->findVariable(R_VARIABLE_PRESSURE);
    if (pressurePos == RConstants::eod)
    {
        pressurePos = this->pModel->addVariable(R_VARIABLE_PRESSURE);
        this->pModel->getVariable(pressurePos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodePressure),
                    RStatistics::findMaximumValue(this->nodePressure));
    }
    RVariable &pressure = this->pModel->getVariable(pressurePos);

    pressure.setApplyType(R_VARIABLE_APPLY_NODE);
    pressure.resize(1,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        pressure.setValue(0,i,this->nodePressure[i]);
    }

    // Velocity
    uint velocityPos = this->pModel->findVariable(R_VARIABLE_VELOCITY);
    if (velocityPos == RConstants::eod)
    {
        velocityPos = this->pModel->addVariable(R_VARIABLE_VELOCITY);
        this->pModel->getVariable(velocityPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumMagnitude(this->elementVelocity),
                    RStatistics::findMaximumMagnitude(this->elementVelocity));
    }
    RVariable &velocity = this->pModel->getVariable(velocityPos);

    velocity.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    velocity.resize(3,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        velocity.setValue(0,i,this->elementVelocity[i][0]);
        velocity.setValue(1,i,this->elementVelocity[i][1]);
        velocity.setValue(2,i,this->elementVelocity[i][2]);
    }

    RLogger::unindent();
}

void RSolverPotential::statistics()
{
    this->printStats(R_VARIABLE_PRESSURE);
    this->printStats(R_VARIABLE_VELOCITY);
    this->processMonitoringPoints();
}

void RSolverPotential::generatePressureNodeBook()
{
    this->nodeBook.resize(this->pModel->getNNodes());
    RBVector disabledPositions(this->nodeBook.size(),false);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        if (!pElementGroup->hasBoundaryCondition(R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT))
        {
            continue;
        }
        for (uint j=0;j<pElementGroup->size();j++)
        {
            const RElement &element = this->pModel->getElement(pElementGroup->get(j));
            for (uint k=0;k<element.size();k++)
            {
                disabledPositions[element.getNodeId(k)] = true;
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
            disabledPositions[i] = true;
        }
    }

    RSolverGeneric::rebuildNodeBook(this->nodeBook,disabledPositions);
}

void RSolverPotential::generateBoundaryVectors(RRVector &elementPressure, RBVector &pressureSetValues)
{
    elementPressure.resize(this->pModel->getNElements());
    elementPressure.fill(0.0);
    pressureSetValues.resize(this->pModel->getNElements());
    pressureSetValues.fill(false);
    this->elementInflowVelocity.resize(this->pModel->getNElements());
    this->elementInflowVelocity.fill(0.0);

    double time = this->pModel->getTimeSolver().getCurrentTime();

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

            if (bc.getType() == R_BOUNDARY_CONDITION_PRESSURE_EXPLICIT)
            {
                uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_PRESSURE);
                if (bcComponentPosition == RConstants::eod)
                {
                    continue;
                }
                double pressure = bc.getComponent(bcComponentPosition).get(time);
                for (uint k=0;k<pElementGroup->size();k++)
                {
                    elementPressure[pElementGroup->get(k)] = pressure;
                    pressureSetValues[pElementGroup->get(k)] = true;
                }
            }
            else if (bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VELOCITY ||
                     bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VOLURATE)
            {
                // Inflow can be prescribed on surfaces only.
                if (pElementGroup->size() == 0 ||
                    !R_ELEMENT_TYPE_IS_SURFACE(this->pModel->getElement(pElementGroup->get(0)).getType()))
                {
                    continue;
                }
                const RSurface *pSurface = static_cast<const RSurface*>(pElementGroup);

                double inflowVelocity = 0.0;
                if (bc.getType() == R_BOUNDARY_CONDITION_INFLOW_VELOCITY)
                {
                    uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_VELOCITY);
                    if (bcComponentPosition == RConstants::eod)
                    {
                        continue;
                    }
                    inflowVelocity = bc.getComponent(bcComponentPosition).get(time);
                }
                else
                {
                    uint bcComponentPosition = bc.findComponentPosition(R_VARIABLE_VOLUME_FLOW_RATE);
                    if (bcComponentPosition == RConstants::eod)
                    {
                        continue;
                    }
                    double area = pSurface->findArea(this->pModel->getNodes(),this->pModel->getElements());
                    if (area < RConstants::eps)
                    {
                        RLogger::warning("Volumetric flow rate boundary condition on entity \'%s\' is ignored - the entity has no area.\n",
                                         pSurface->getName().toUtf8().constData());
                        continue;
                    }
                    inflowVelocity = bc.getComponent(bcComponentPosition).get(time) / area;
                }

                for (uint k=0;k<pSurface->size();k++)
                {
                    this->elementInflowVelocity[pSurface->get(k)] = inflowVelocity;
                }
            }
        }
    }
}

RR3Vector RSolverPotential::getElementGravity(uint elementID) const
{
    return RR3Vector(this->elementGravityX[elementID],
                     this->elementGravityY[elementID],
                     this->elementGravityZ[elementID]);
}

void RSolverPotential::assemblyMatrix(unsigned int elementID, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    const RElement &element = this->pModel->getElement(elementID);

    RRMatrix Ae(Ke);
    RRVector be(fe);

    // Apply explicit boundary conditions.
    for (uint m=0;m<element.size();m++)
    {
        uint position;
        uint nodeID = element.getNodeId(m);
        if (!this->nodeBook.getValue(nodeID,position))
        {
            for (uint n=0;n<element.size();n++)
            {
                be[n] -= Ae[n][m] * this->nodePressure[nodeID];
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

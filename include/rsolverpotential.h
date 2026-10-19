#ifndef RSOLVERPOTENTIAL_H
#define RSOLVERPOTENTIAL_H

#include "rsolvergeneric.h"

//! Steady-state flow of an incompressible fluid through porous media.
//! The fluid velocity follows Darcy's law v = -(k/mu)*(grad(p) - rho*g),
//! which together with continuity div(v) = 0 yields a potential (Laplace type)
//! equation for the pressure p.
class RSolverPotential : public RSolverGeneric
{

    protected:

        //! Node pressure.
        RRVector nodePressure;
        //! Element Darcy (superficial) velocity.
        std::vector<RR3Vector> elementVelocity;
        //! Element permeability.
        RRVector elementPermeability;
        //! Element dynamic viscosity.
        RRVector elementViscosity;
        //! Element density.
        RRVector elementDensity;
        //! Element hydraulic conductivity k/mu.
        RRVector elementConductivity;
        //! Conductivity the matrix system is normalized by.
        double conductivityScale;
        //! Element gravitational acceleration.
        RRVector elementGravityX;
        RRVector elementGravityY;
        RRVector elementGravityZ;
        //! Element inflow velocity (normal, positive into the domain).
        RRVector elementInflowVelocity;

    public:

        //! Constructor.
        explicit RSolverPotential(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverPotential();

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        //! Initialize solver.
        void initialize() override;

        //! Update scales.
        void updateScales() override;

        //! Recover previously computed results.
        void recover() override;

        //! Prepare solver.
        void prepare() override;

        //! Run matrix solver.
        void solve() override;

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Generate node book - only explicit pressure fixes node values.
        //! Inflow conditions prescribe a flux and therefore stay free.
        void generatePressureNodeBook();

        //! Generate element vectors of explicit pressure and inflow velocity.
        void generateBoundaryVectors(RRVector &elementPressure, RBVector &pressureSetValues);

        //! Return element gravitational acceleration.
        RR3Vector getElementGravity(uint elementID) const;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

};

#endif // RSOLVERPOTENTIAL_H

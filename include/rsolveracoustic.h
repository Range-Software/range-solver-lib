#ifndef RSOLVERACOUSTIC_H
#define RSOLVERACOUSTIC_H

#include "rsolvergeneric.h"

class RSolverAcoustic : public RSolverGeneric
{

    protected:

        //! Element modulus of elasticity vector.
        RRVector elementElasticityModulus;
        //! Element density vector.
        RRVector elementDensity;
        //! Element damping factor.
        RRVector elementDampingFactor;

        //! Node velocity potential.
        RRVector nodeVelocityPotential;
        //! Node velocity potential - old.
        RRVector nodeVelocityPotentialOld;
        //! Node velocity potential - velocity.
        RRVector nodeVelocityPotentialVelocity;
        //! Node velocity potential - acceleration.
        RRVector nodeVelocityPotentialAcceleration;
        //! Node acoustic pressure.
        RRVector nodeAcousticPressure;
        //! Element acoustic pressure.
        RSolverCartesianVector<RRVector> elementAcousticParticleVelocity;

    public:

        //! Constructor.
        explicit RSolverAcoustic(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverAcoustic();

        //! Check if solver has converged.
        bool hasConverged() const override;

    protected:

        double findSoundSpeedScale() const;

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

        //! Process absorbing boundary.
        void processAbsorbingBoundary();

        //! Process acoustic pressure.
        void processAcousticPressure();

        //! Process acoustic particle velocity.
        void processAcousticParticleVelocity();

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ce, const RRMatrix &Ke, const RRVector &fe);

        //! Find absorbing boundary nodes.
        std::vector<bool> findAbsorbingBoundaryNodes() const;

        //! Find absorbing boundary elements.
        std::vector<bool> findAbsorbingBoundaryElements() const;

        //! Find absorbing boundary nodes.
        std::vector<RR3Vector> findAbsorbingBoundaryNormals(const std::vector<bool> &absorbingBoundaryNodes) const;

};

#endif // RSOLVERACOUSTIC_H

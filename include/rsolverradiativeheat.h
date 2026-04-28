#ifndef RSOLVERRADIATIVEHEAT_H
#define RSOLVERRADIATIVEHEAT_H

#include <rbl_rvector.h>

#include "rsolvergeneric.h"


class RSolverRadiativeHeat : public RSolverGeneric
{

    protected:

        //! Element temperature.
        RRVector elementTemperature;
        //! Element heat vector.
        RRVector elementRadiativeHeat;
        //! Patch heat vector.
        RRVector patchHeat;
        //! View-factor matrix.
        RViewFactorMatrix viewFactorMatrix;
        //! Patch heat norm.
        double patchHeatNorm;
        //! Old patch heat norm.
        double oldPatchHeatNorm;

    public:

        //! Constructor.
        explicit RSolverRadiativeHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverRadiativeHeat();

        //! Check if solver has converged.
        bool hasConverged() const;

    protected:

        //! Find temperature scale.
        double findTemperatureScale() const;

        //! Initialize solver.
        void initialize() override;

        //! Update scales.
        void updateScales();

        //! Recover previously computed results.
        void recover();

        //! Prepare view-factors.
        void prepareViewFactors();

        //! Prepare solver.
        void prepare();

        //! Run matrix solver.
        void solve();

        //! Process solver results.
        void process();

        //! Store solver results.
        void store();

        //! Process statistics.
        void statistics();

        //! Check if view-factor header correspond with input.
        bool checkViewFactorHeader(const RViewFactorMatrixHeader &viewFactorMatrixHeader) const;

};

#endif // RSOLVERRADIATIVEHEAT_H

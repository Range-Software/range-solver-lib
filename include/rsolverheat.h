#ifndef RSOLVERHEAT_H
#define RSOLVERHEAT_H

#include <rbl_r3vector.h>

#include "rsolvergeneric.h"

class RSolverHeat : public RSolverGeneric
{

    protected:

        //! Element heat capacity vector.
        RRVector elementCapacity;
        //! Element thermal conduction vector.
        RRVector elementConduction;
        //! Element density vector.
        RRVector elementDensity;
        //! Node temperature.
        RRVector nodeTemperature;
        //! Element heat vector.
        RRVector elementHeat;
        //! Element radiation heat vector.
        RRVector elementRadiativeHeat;
        //! Element joule heat.
        RRVector elementJouleHeat;
        //! Element heat flux vector.
        std::vector<RR3Vector> elementHeatFlux;

    public:

        //! Constructor.
        explicit RSolverHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverHeat();

        //! Check if solver has converged.
        bool hasConverged() const;

    protected:

        //! Find temperature scale.
        double findTemperatureScale() const;

        //! Update scales.
        void updateScales();

        //! Recover previously computed results.
        void recover();

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

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe);

        //! Get simple convection BC values.
        bool getSimpleConvection(const RElementGroup &elementGroup, double &htc, double &htt);

        //! Get forced convection BC values.
        bool getForcedConvection(const RElementGroup &elementGroup, double &htc, double &htt);

        //! Get natural convection BC values.
        bool getNaturalConvection(const RElementGroup &elementGroup, unsigned int elementId, double &htc, double &htt);

};

#endif // RSOLVERHEAT_H

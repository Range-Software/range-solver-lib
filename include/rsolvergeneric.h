#ifndef RSOLVERGENERIC_H
#define RSOLVERGENERIC_H

#include <QString>

#include <vector>

#include <rbl_book.h>
#include <rml_sparse_matrix.h>

#include "rlocalrotation.h"
#include "rscales.h"
#include "rsolvershareddata.h"

template <typename T>
struct RSolverCartesianVector
{
    public:

        T x, y, z;

};

class RSolverGeneric
{

    protected:

        //! Vacuum permittivity.
        static const double e0;
        //! Stefan–Boltzmann constant.
        static const double sigma;

    protected:

        //! Indicator whether the mesh has changed.
        bool meshChanged;

        //! Problem type.
        RProblemType problemType;
        //! Pointer to model object.
        RModel *pModel;
        //! Model file name.
        QString modelFileName;
        //! Convergence file.
        QString convergenceFileName;
        //! Matrix M (modal analysis).
        RSparseMatrix M;
        //! Matrix A.
        RSparseMatrix A;
        //! Vector x.
        RRVector x;
        //! Vector b.
        RRVector b;
        //! Node book.
        RBook nodeBook;
        //! Local rotations.
        std::vector<RLocalRotation> localRotations;
        //! Element temperature vector.
        RRVector elementTemperature;
        //! Scales.
        RScales scales;
        //! Pointer to shared data container.
        RSolverSharedData *pSharedData;
        //! First run indicator.
        bool firstRun;
        //! Task iteration counter.
        uint taskIteration;
        //! Computable elements array.
        RBVector computableElements;
        //! Includable elements array (not included in coefficient matrix).
        RBVector includableElements;
        //! Surface element inward orientation (if normal is pointing inside computable volume element).
        RBVector inwardElements;

    public:

        //! Constructor.
        RSolverGeneric(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Copy constructor (deleted - solvers should not be copied).
        RSolverGeneric(const RSolverGeneric &genericSolver) = delete;

        //! Assignment operator (deleted - solvers should not be copied).
        RSolverGeneric & operator =(const RSolverGeneric &genericSolver) = delete;

        //! Destructor.
        virtual ~RSolverGeneric();

        //! Run solver.
        void run(bool firstRun, uint taskIteration);

        //! Get mesh changed.
        bool getMeshChanged() const;

        //! Set mesh changed.
        void setMeshChanged(bool meshChanged);

        //! Check if solver has converged.
        virtual bool hasConverged() const = 0;

        //! Update old records.
        static void updateOldRecords(const RTimeSolver &rTimeSolver, const QString &modelFileName);

    protected:

        //! Find element sizes.
        RRVector findElementSizes() const;

        //! Find element scale.
        double findElementScale(bool onlyVolumes = false) const;

        //! Find mesh scale.
        double findMeshScale() const;

        //! Update local rotations.
        void updateLocalRotations();

        //! Clear shared data.
        virtual void clearSharedData();

        //! Store shared data.
        virtual void storeSharedData();

        //! Recover shared data.
        virtual void recoverSharedData();

        //! Update scales.
        virtual void updateScales() = 0;

        //! Recover previously computed results.
        virtual void recover() = 0;

        //! Prepare solver.
        virtual void prepare() = 0;

        //! Run matrix solver.
        virtual void solve() = 0;

        //! Process solver results.
        virtual void process() = 0;

        //! Store solver results.
        virtual void store() = 0;

        //! Process statistics.
        virtual void statistics() = 0;

        //! Write results.
        void writeResults();

        //! Apply displacement if possible.
        void applyDisplacement();

        //! Remove displacement if possible.
        void removeDisplacement();

        //! Generate node book.
        void generateNodeBook(RProblemType problemType);

        //! Rebuild node book from disabled positions.
        static void rebuildNodeBook(RBook &nodeBook, const RBVector &disabledPositions);

        //! Generate material element vector.
        void generateMaterialVecor(RMaterialProperty::Type materialPropertyType, RRVector &materialPropertyValues) const;

        //! Generate variable element vector from applied environment/initial/boundary condition.
        void generateVariableVector(RVariableType variableType,
                                    RRVector &variableValues,
                                    RBVector &setValues,
                                    bool boundaryConditions,
                                    bool initialConditions,
                                    bool environmentConditions,
                                    bool onlyExplicitBcs = false) const;

        //! Recover variable.
        void recoverVariable(RVariableType variableType,
                             RVariableApplyType applyType,
                             uint expectedSize,
                             uint vectorPosition,
                             RRVector &variableValues,
                             double defaultValue) const;

        //! Sunc values with shared data.
        void syncShared(const QString &keyName, RRVector &values);

        //! Find computable elements based on assigned material properties or has boundary condition assigned.
        void findComputableElements(RProblemType problemType);

        //! Find includable elements (not computed in coefficient matrix) - all points are computable.
        void findIncludableElements();

        //! Find inward surface elements (normal is pointing inside computable element).
        void findInwardElements();

        //! Process monitoring points.
        void processMonitoringPoints() const;

        //! Print results statistics.
        void printStats(RVariableType variableType) const;

};

#endif // RSOLVERGENERIC_H

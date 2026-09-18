#ifndef RSOLVERMAGNETOSTATICS_H
#define RSOLVERMAGNETOSTATICS_H

#include <array>
#include <vector>

#include "rsolvergeneric.h"

//! Magneto-static solver.
//!
//! The magnetic flux density is evaluated directly from the current density
//! computed by the electro-static task with the Biot-Savart law. There is no
//! system of equations to solve and no boundary condition to impose - the
//! condition that the field vanishes at infinity is built into the integral.
class RSolverMagnetostatics : public RSolverGeneric
{

    public:

        //! Current carrying source cell of the Biot-Savart integral - a
        //! segment, a triangle or a tetrahedron with a uniform current.
        struct Source
        {
            //! Number of vertices - 2 for a segment, 3 for a triangle, 4 for a tetrahedron.
            uint nVertices;
            //! Vertex coordinates.
            std::array<std::array<double,3>,4> x;
            //! Current carried by the cell:
            //! segment     - current J * cross area [A], only its component
            //!               along the segment is used,
            //! triangle    - sheet current J * thickness [A/m],
            //! tetrahedron - current density J [A/m^2].
            std::array<double,3> j;
            //! Centroid.
            std::array<double,3> center;
            //! Square of the longest edge.
            double size2;
            //! Length, area or volume.
            double measure;

            //! Create a source cell and compute its centroid, size and measure.
            static Source create(uint nVertices, const std::array<std::array<double,3>,4> &x, const std::array<double,3> &j);
        };

    protected:

        //! Element current density.
        RSolverCartesianVector<RRVector> elementCurrentDensity;
        //! Current carrying source cells.
        std::vector<Source> sources;
        //! Node magnetic field.
        RSolverCartesianVector<RRVector> nodeMagneticField;

    public:

        //! Constructor.
        explicit RSolverMagnetostatics(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverMagnetostatics();

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Magnetic flux density [T] induced at point p by a single source cell.
        //! Cells close to p are integrated exactly, distant ones with the
        //! midpoint rule.
        static std::array<double,3> findSourceField(const Source &source, const std::array<double,3> &p);

        //! Magnetic flux density [T] induced at point p by a single source cell,
        //! always integrated exactly.
        static std::array<double,3> findSourceFieldExact(const Source &source, const std::array<double,3> &p);

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

};

#endif // RSOLVERMAGNETOSTATICS_H

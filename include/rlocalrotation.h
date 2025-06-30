#ifndef RLOCALROTATION_H
#define RLOCALROTATION_H

#include <rbl_rmatrix.h>

class RLocalRotation
{

    protected:

        //! Rotation is active.
        bool active;
        //! Rotation matrix.
        RRMatrix R;
        //! Inverse rotation matrix.
        RRMatrix iR;

    private:

        //! Internal initialization function.
        void _init(const RLocalRotation *pLocalRotation = nullptr);

    public:

        //! Constructor.
        RLocalRotation();

        //! Copy constructor.
        RLocalRotation(const RLocalRotation &localRotation);

        //! Destructor.
        ~RLocalRotation();

        //! Assignment operator.
        RLocalRotation &operator =(const RLocalRotation &localRotation);

        //! Return true if local rotation is active.
        bool isActive(void) const;

        //! Activate local rotation.
        void activate(const RRMatrix &R);

        //! Deactivate local rotation.
        void deactivate(void);

        //! Return const reference to rotation matrix.
        inline const RRMatrix &getR(void) const
        {
            return this->R;
        }

        //! Return const reference to inverse rotation matrix.
        inline const RRMatrix &getInverseR(void) const
        {
            return this->iR;
        }

        //! Rotate vector.
        void rotateResultsVector(RRVector &v) const;

};

#endif // RLOCALROTATION_H

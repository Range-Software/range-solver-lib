#ifndef RHEMICUBEPIXEL_H
#define RHEMICUBEPIXEL_H

#include <rbl_r3vector.h>
#include <rbl_utils.h>

class RHemiCubePixel
{

    protected:

        //! Pixel position.
        RR3Vector position;
        //! Pixel color (0 = no color).
        uint color;
        //! Pixel depth.
        double depth;
        //! Weight factor.
        double weight;

    private:

        //! Internal initialization function.
        void _init(const RHemiCubePixel *pHemiCubePixel = nullptr);

    public:

        //! Constructor.
        RHemiCubePixel();

        //! Copy constructor.
        RHemiCubePixel(const RHemiCubePixel &hemiCubePixel);

        //! Destructor.
        ~RHemiCubePixel();

        //! Assignment operator.
        RHemiCubePixel &operator =(const RHemiCubePixel &hemiCubePixel);

        //! Return const reference to pixel position.
        const RR3Vector &getPosition() const;

        //! Return reference to pixel position.
        RR3Vector &getPosition();

        //! Set pixel position.
        void setPosition(const RR3Vector &position);

        //! Return pixel color.
        uint getColor() const;

        //! Set pixel color.
        void setColor(uint color);

        //! Return pixel depth.
        double getDepth() const;

        //! Set pixel depth.
        void setDepth(double depth);

        //! Return pixel weight.
        double getWeight() const;

        //! Set pixel weight.
        void setWeight(double weight);

};

#endif // RHEMICUBEPIXEL_H

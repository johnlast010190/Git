/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2010-2012 Engys Ltd.

 @author Jacques Papper
 \*---------------------------------------------------------------------------*/

#ifndef _CURVEFITTER2D_H
#define    _CURVEFITTER2D_H


#include <vector>
#include "interpolation/curveFitting/pieceWiseCubicCurveFitter.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "fields/Fields/vector2DField/vector2DFieldFwd.H"
#include "containers/Lists/List/List.H"
#include "primitives/Tuple2/Tuple2.H"

namespace Foam
{

    class curveFitter2D
    {
    public:

        curveFitter2D();
        curveFitter2D(const List< Tuple2< scalar, vector2DField >>  & allCurves);
        curveFitter2D(const curveFitter2D& orig);
        virtual ~curveFitter2D();

        scalar interpolate(const scalar & input1, const scalar & input2);
        void addCurve(const scalar & X, const vector2DField & curve);
        void addAllCurves(const List< Tuple2< scalar, vector2DField >>  & allCurves);

    protected:

        std::vector< pieceWiseCubicCurveFitter > alldata_;
        std::vector< scalar > allDataX_;


    };
};
#endif    /* _CURVEFITTER2D_H */

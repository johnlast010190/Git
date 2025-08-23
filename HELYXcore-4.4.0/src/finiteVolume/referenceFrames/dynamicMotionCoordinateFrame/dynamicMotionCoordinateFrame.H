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
    (c) 2023 Engys Ltd.

Class
    Foam::dynamicMotionCoordinateFrame

Description
    Moving frame of reference based on motion function.

Usage
    Example of the vehicle frame specification:
    \verbatim
    <referenceFrameName>
    {
        type            dynamicMotionCoordinateFrame;
        motionFunction
        {
            type vehicleMotion;
            vehicleMotionCoeffs
            {
                vehicleSpeed    constant 50;
                invTurnRadius   table ( ( 0 0.02 ) ( 2000 0.02 ) );
            }
        }
        coordinateSystem
        {
            type            cartesian;
            origin          ( 0 50 -0.05 );
            e1              ( -1 0 0 );
            e2              ( 0 -1 0 );
        }
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionCoordinateFrame_H
#define dynamicMotionCoordinateFrame_H

#include "referenceFrames/motionCoordinateFrame/motionCoordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dynamicMotionCoordinateFrame Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionCoordinateFrame
:
    public motionCoordinateFrame
{

public:

    //- Runtime type information
    TypeName("dynamicMotionCoordinateFrame");

    //- Constructors
        dynamicMotionCoordinateFrame
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& frameName
        );

    //- Destructor
        virtual ~dynamicMotionCoordinateFrame();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//end namespace FOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

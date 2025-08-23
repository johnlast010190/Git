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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::solidBodyMotionFunctions::crossWind


SourceFiles
    crossWind.C

\*---------------------------------------------------------------------------*/

#ifndef crossWind_H
#define crossWind_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class crossWind Declaration
\*---------------------------------------------------------------------------*/

class crossWind
:
    public solidBodyMotionFunction
{
    // Private data

        //- Linear velocity
        autoPtr<Function1<vector>> velocity_;

        //- Name of the frame for rotation
        const word frameName_;

        //- Frame for rotation
        const coordinateFrame& frameForRotation_;


public:

    //- Runtime type information
    TypeName("crossWind");


    // Constructors

        //- Construct from components
        crossWind
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Disallow default bitwise copy construction
        crossWind(const crossWind&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new crossWind(time_, SBMFCoeffs_, frameName_)
            );
        }


    //- Destructor
    virtual ~crossWind();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t1,
            const scalar t2
        ) const override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;

        //- Return the motion velocity
        virtual vectorTuple acceleration() const override
        {
            return vectorTuple(Zero, Zero);
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const crossWind&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

Class
    Foam::solidBodyMotionFunctions::tabulated6DoFMotion

Description
    Tabulated 6DoF motion function.

    Obtained by interpolating tabulated data for surge (x-translation),
    sway (y-translation), heave (z-translation), roll (rotation about x),
    pitch (rotation about y) and yaw (rotation about z).

SourceFiles
    tabulated6DoFMotion.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated6DoFMotion_H
#define tabulated6DoFMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "primitives/Vector2D/vector2DVector/vector2DVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                    Class tabulated6DoFMotion Declaration
\*---------------------------------------------------------------------------*/

class tabulated6DoFMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Time data file name read from dictionary
        fileName timeDataFileName_;

        //- Field of times
        scalarField times_;

        //- Field of translation and rotation "vectors"
        Field<vector2DVector> values_;

        //- Center of gravity read from dictionary
        mutable vector origin_;


public:

    //- Runtime type information
    TypeName("tabulated6DoFMotion");


    // Constructors

        //- Construct from components
        tabulated6DoFMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        tabulated6DoFMotion
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Disallow default bitwise copy construction
        tabulated6DoFMotion(const tabulated6DoFMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new tabulated6DoFMotion(SBMFCoeffs_, time_, frameName_)
            );
        }


    //- Destructor
    virtual ~tabulated6DoFMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t1,
            const scalar t2
        ) const override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tabulated6DoFMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

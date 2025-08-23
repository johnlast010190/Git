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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

Class
    Foam::solidBodyMotionFunctions::rotatingStepMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    rotatingStepMotion.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingStepMotion_H
#define rotatingStepMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingStepMotion Declaration
\*---------------------------------------------------------------------------*/

class rotatingStepMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Origin of the axis
        mutable vector origin_;

        //- Axis vector
        mutable vector axis_;

        //- Step angle (degrees)
        scalar theta_;

        //- Initial angle at index 0
        mutable scalar theta0_;

        //- Step period
        scalar stepPeriod_;

        //- Previous step update
        mutable label nStepsDone_;


public:

    //- Runtime type information
    TypeName("rotatingStepMotion");


    // Constructors

        //- Construct from components
        rotatingStepMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        rotatingStepMotion
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Disallow copy construct
        rotatingStepMotion(const rotatingStepMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new rotatingStepMotion(SBMFCoeffs_, time_, frameName_)
            );
        }

    //- Destructor
    virtual ~rotatingStepMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rotatingStepMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

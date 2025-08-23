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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::solidBodyMotionFunctions::vehicleMotion

Description
    Motion should now support use in both as just a momentum source and for
    dynamic mesh cases in incremental mode.

SourceFiles
    vehicleMotion.C

\*---------------------------------------------------------------------------*/

#ifndef vehicleMotion_H
#define vehicleMotion_H

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
                          Class vehicleMotion Declaration
\*---------------------------------------------------------------------------*/

class vehicleMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of rotation
        mutable vector CofR_;

        //- Centre of rotation
        mutable vector CofR0_;

        //- vehicle linear speed
        autoPtr<Function1<scalar>> vehicleSpeed_;

        //- Turn centroid in vehicle frame
        autoPtr<Function1<scalar>> turnRadius_;

        //- Inverse turn centroid in vehicle frame
        autoPtr<Function1<scalar>> invTurnRadius_;

        //- Vehicle angular speed
        mutable vector omega_;

        //- Linear velocity
        mutable vector linearVelocity_;


public:

    //- Runtime type information
    TypeName("vehicleMotion");


    // Constructors

        //- Construct from components
        vehicleMotion
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Disallow default bitwise copy construction
        vehicleMotion(const vehicleMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new vehicleMotion(time_, SBMFCoeffs_, frameName_)
            );
        }


    //- Destructor
    virtual ~vehicleMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const override;

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t1,
            const scalar t2
        ) const override
        {
            NotImplemented;
        }

        //- Same update but allows even from the const functions
        void updateMotionCoeffs() const;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;

        //- Return centre of rotation
        virtual const vector& CofR() const override
        {
            return CofR_;
        }

        //- Return centre of rotation
        virtual const vector& CofR0() const override
        {
            return CofR0_;
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const vehicleMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

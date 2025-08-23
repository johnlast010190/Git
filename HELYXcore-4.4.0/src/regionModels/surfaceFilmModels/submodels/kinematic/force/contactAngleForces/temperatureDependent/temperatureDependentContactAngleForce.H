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
    (c) 2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::temperatureDependentContactAngleForce

Description
    Temperature dependent contact angle force

    The contact angle in degrees is specified as a Foam::Function1 type, to
    enable the use of, e.g.  contant, polynomial, table values.

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForce
    Foam::Function1Types

SourceFiles
    temperatureDependentContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentContactAngleForce_H
#define temperatureDependentContactAngleForce_H

#include "surfaceFilmModels/submodels/kinematic/force/contactAngleForces/contactAngleForce/contactAngleForce.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
            Class temperatureDependentContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Contact angle function
        autoPtr<Function1<scalar>> thetaPtr_;


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("temperatureDependentContactAngle");


    // Constructors

        //- Construct from surface film model
        temperatureDependentContactAngleForce
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        temperatureDependentContactAngleForce
        (
            const temperatureDependentContactAngleForce&
        ) = delete;


    //- Destructor
    virtual ~temperatureDependentContactAngleForce();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const temperatureDependentContactAngleForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

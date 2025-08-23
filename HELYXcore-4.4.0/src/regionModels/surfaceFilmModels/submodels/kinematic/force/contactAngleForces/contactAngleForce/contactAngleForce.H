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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::contactAngleForce

Description
    Base-class for film contact angle force models.

    The effect of the contact angle force can be ignored over a specified
    distance from patches.

SourceFiles
    contactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef contactAngleForce_H
#define contactAngleForce_H

#include "surfaceFilmModels/submodels/kinematic/force/force/force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class contactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class contactAngleForce
:
    public force
{
    // Private Data

        //- Coefficient applied to the contact angle force
        scalar Ccf_;

        //- Mask over which force is applied
        volScalarField mask_;


    // Private member functions

        //- Initialise
        void initialise();


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const = 0;


public:

    //- Runtime type information
    TypeName("contactAngle");


    // Constructors

        //- Construct from surface film model
        contactAngleForce
        (
            const word& typeName,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        contactAngleForce(const contactAngleForce&) = delete;


    //- Destructor
    virtual ~contactAngleForce();


    // Member Functions

        //- Correct
        virtual tmp<fvVectorMatrix> correct(volVectorField& U);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const contactAngleForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

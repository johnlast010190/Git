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
    (c) 2017-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::distributionContactAngleForce

Description
    PDF distribution based film contact angle force

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForce
    Foam::distributionModel

SourceFiles
    distributionContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef distributionContactAngleForce_H
#define distributionContactAngleForce_H

#include "surfaceFilmModels/submodels/kinematic/force/contactAngleForces/contactAngleForce/contactAngleForce.H"
#include "distributionModel/distributionModel.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                Class distributionContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class distributionContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;

        //- Cached contact angle field
        mutable volScalarField theta_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("distributionContactAngle");


    // Constructors

        //- Construct from surface film model
        distributionContactAngleForce
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        distributionContactAngleForce
        (
            const distributionContactAngleForce&
        ) = delete;


    //- Destructor
    virtual ~distributionContactAngleForce();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const distributionContactAngleForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::standardRadiation

Description
    Standard radiation model

SourceFiles
    standardRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef standardRadiation_H
#define standardRadiation_H

#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class standardRadiation Declaration
\*---------------------------------------------------------------------------*/

class standardRadiation
:
    public radiationModel
{
    // Private data

        //- Radiative incident flux mapped from the primary region [kg/s^3]
        volScalarField::Internal qinFilm_;

        //- Beta coefficient
        scalar beta_;

        //- Bar(kappa) coefficient
        scalar kappaBar_;


public:

    //- Runtime type information
    TypeName("standardRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        standardRadiation
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        standardRadiation(const standardRadiation&) = delete;


    //- Destructor
    virtual ~standardRadiation();


    // Member Functions

        //- Correct
        virtual void correct();

        //- Return the radiation sensible enthalpy source
        //  Also updates qrNet
        virtual tmp<volScalarField::Internal> Shs();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const standardRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

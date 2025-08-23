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
    (c) 2012-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/primaryRadiation/primaryRadiation.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primaryRadiation, 0);
addToRunTimeSelectionTable(radiationModel, primaryRadiation, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primaryRadiation::primaryRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    radiationModel(typeName, film, dict),
    qinFilm_
    (
        IOobject
        (
            "qin",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

primaryRadiation::~primaryRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void primaryRadiation::correct()
{
    const volScalarField& qinPrimary
    (
        film().primaryMesh().lookupObject<volScalarField>("qin")
    );

    // Map the primary-side radiative flux to the region internal field
    film().toRegion(qinFilm_, qinPrimary.boundaryField());
}


tmp<volScalarField::Internal> primaryRadiation::Shs()
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Shs", typeName),
        qinFilm_*filmModel_.coverage()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

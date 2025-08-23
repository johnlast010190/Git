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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/constantRadiation/constantRadiation.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantRadiation, 0);
addToRunTimeSelectionTable(radiationModel, constantRadiation, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantRadiation::constantRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    radiationModel(typeName, film, dict),
    qrConst_
    (
        IOobject
        (
            IOobject::modelName("qrConst", typeName),
            film.time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    ),
    mask_
    (
        IOobject
        (
            IOobject::modelName("mask", typeName),
            film.time().timeName(),
            film.regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimless, 1)
    ),
    absorptivity_(coeffDict_.lookup<scalar>("absorptivity")),
    timeStart_(coeffDict_.lookup<scalar>("timeStart")),
    duration_(coeffDict_.lookup<scalar>("duration"))
{
    mask_ = pos0(mask_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantRadiation::~constantRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void constantRadiation::correct()
{}


tmp<volScalarField::Internal> constantRadiation::Shs()
{
    const scalar time = film().time().value();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        return volScalarField::Internal::New
        (
            IOobject::modelName("Shs", typeName),
            mask_*qrConst_*filmModel_.coverage()*absorptivity_
        );
    }
    else
    {
        return volScalarField::Internal::New
        (
            IOobject::modelName("Shs", typeName),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), 0)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

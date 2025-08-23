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
    (c) 2011-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/noFilm/noFilm.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noFilm, 0);
addToRunTimeSelectionTable(surfaceFilmModel, noFilm, mesh);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noFilm::noFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType
)
:
    surfaceFilmModel(),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noFilm::~noFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar noFilm::CourantNumber() const
{
    return 0;
}


tmp<volScalarField::Internal> noFilm::Srho() const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Srho", typeName),
        mesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    );
}


tmp<volScalarField::Internal> noFilm::Srho(const label i) const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Srho(" + Foam::name(i) + ")", typeName),
        mesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    );
}


tmp<volVectorField::Internal> noFilm::SU() const
{
    return volVectorField::Internal::New
    (
        IOobject::modelName("SU", typeName),
        mesh_,
        dimensionedVector(dimMass/dimVolume/dimTime, Zero)
    );
}


tmp<volScalarField::Internal> noFilm::Sh() const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Sh", typeName),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


void noFilm::evolve()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

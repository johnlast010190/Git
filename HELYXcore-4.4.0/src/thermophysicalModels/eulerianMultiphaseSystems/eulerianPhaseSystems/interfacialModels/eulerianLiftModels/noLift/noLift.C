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
    (c) 2014-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "noLift.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianLiftModels
{
    defineTypeNameAndDebug(noLift, 0);
    addToRunTimeSelectionTable(eulerianLiftModel, noLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::noLift::noLift
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianLiftModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::noLift::~noLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianLiftModels::noLift::Cl() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();
    return volScalarField::New("Cl", mesh, dimensionedScalar(dimless, 0));
}


Foam::tmp<Foam::volVectorField> Foam::eulerianLiftModels::noLift::F() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();
    return volVectorField::New
    (
        "noLift:F",
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::eulerianLiftModels::noLift::Ff() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();
    return surfaceScalarField::New
    (
        "noLift:Ff",
        mesh,
        dimensionedScalar(dimF*dimArea, 0)
    );
}


// ************************************************************************* //

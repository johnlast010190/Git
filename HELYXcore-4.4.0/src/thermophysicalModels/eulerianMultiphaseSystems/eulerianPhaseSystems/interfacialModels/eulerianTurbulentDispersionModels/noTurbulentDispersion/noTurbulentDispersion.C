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

#include "noTurbulentDispersion.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianTurbulentDispersionModels
{
    defineTypeNameAndDebug(noTurbulentDispersion, 0);
    addToRunTimeSelectionTable
    (
        eulerianTurbulentDispersionModel,
        noTurbulentDispersion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::noTurbulentDispersion::noTurbulentDispersion
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianTurbulentDispersionModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::noTurbulentDispersion::
~noTurbulentDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianTurbulentDispersionModels::noTurbulentDispersion::D() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();
    return volScalarField::New("zero", mesh, dimensionedScalar(dimD, 0));
}


Foam::tmp<Foam::volVectorField>
Foam::eulerianTurbulentDispersionModels::noTurbulentDispersion::F() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();
    return volVectorField::New("zero", mesh, dimensionedVector(dimF, Zero));
}


// ************************************************************************* //

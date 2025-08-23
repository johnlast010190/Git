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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "noWallDamping.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianWallDampingModels
{
    defineTypeNameAndDebug(noWallDamping, 0);
    addToRunTimeSelectionTable
    (
        eulerianWallDampingModel,
        noWallDamping,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianWallDampingModels::noWallDamping::noWallDamping
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianWallDampingModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianWallDampingModels::noWallDamping::~noWallDamping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianWallDampingModels::noWallDamping::damp
(
    const tmp<volScalarField>& Cl
) const
{
    return Cl;
}


Foam::tmp<Foam::volVectorField>
Foam::eulerianWallDampingModels::noWallDamping::damp
(
    const tmp<volVectorField>& F
) const
{
    return F;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::eulerianWallDampingModels::noWallDamping::damp
(
    const tmp<surfaceScalarField>& Ff
) const
{
    return Ff;
}


// ************************************************************************* //

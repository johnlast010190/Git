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

#include "noWallLubrication.H"
#include "phasePair/phasePair/phasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(noWallLubrication, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        noWallLubrication,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::noWallLubrication::noWallLubrication
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallLubricationModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::noWallLubrication::~noWallLubrication()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::noWallLubrication::Fi() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return volVectorField::New
    (
        "noWallLubrication-Fi",
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::noWallLubrication::F() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return volVectorField::New
    (
        "noWallLubrication-F",
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


// ************************************************************************* //

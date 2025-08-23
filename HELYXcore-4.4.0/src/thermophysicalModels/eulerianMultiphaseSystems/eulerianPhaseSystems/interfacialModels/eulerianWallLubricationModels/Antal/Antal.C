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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Antal.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianWallLubricationModels
{
    defineTypeNameAndDebug(Antal, 0);
    addToRunTimeSelectionTable
    (
        eulerianWallLubricationModel,
        Antal,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModels::Antal::Antal
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianWallLubricationModel(dict, pair),
    Cw1_("Cw1", dimless, dict),
    Cw2_("Cw2", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModels::Antal::~Antal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::eulerianWallLubricationModels::Antal::Fi() const
{
    volVectorField Ur(pair_.Ur());

    const volVectorField& n(nWall());

    return zeroGradWalls
    (
        max
        (
            dimensionedScalar(dimless/dimLength, 0),
            Cw1_/pair_.dispersed().d() + Cw2_/yWall()
        )
       *pair_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //

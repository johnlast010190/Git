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

#include "sphericalHeatTransfer.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianHeatTransferModels
{
    defineTypeNameAndDebug(sphericalHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        eulerianHeatTransferModel,
        sphericalHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianHeatTransferModels::sphericalHeatTransfer::sphericalHeatTransfer
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianHeatTransferModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianHeatTransferModels::sphericalHeatTransfer::~sphericalHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianHeatTransferModels::sphericalHeatTransfer::K
(
    const scalar residualAlpha
) const
{
    return
        60.0
       *max(pair_.dispersed().volFrac(), residualAlpha)
       *pair_.continuous().kappa()
       /sqr(pair_.dispersed().d());
}


// ************************************************************************* //

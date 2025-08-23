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
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "noMassTransfer.H"
#include "solverObjects/disperseEulerian/phase/phase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{
namespace massTransferModels
{
    defineTypeNameAndDebug(noMassTransfer, 0);
    addToRunTimeSelectionTable(massTransferModel, noMassTransfer, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decoupledEulerian::massTransferModels::noMassTransfer::noMassTransfer
(
    const dictionary& dict,
    const phase& phase
)
:
    massTransferModel(dict, phase)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decoupledEulerian::massTransferModels::noMassTransfer::~noMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::decoupledEulerian::massTransferModels::noMassTransfer::mDot() const
{
    return volScalarField::New
    (
        "noMassTransfer-mDot",
        phase_.alphad().mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::decoupledEulerian::massTransferModels::noMassTransfer::Sh() const
{
    return volScalarField::New
    (
        "noMassTransfer-Sh",
        phase_.alphad().mesh(),
        dimensionedScalar(dimless, Zero)
    );
}

// ************************************************************************* //

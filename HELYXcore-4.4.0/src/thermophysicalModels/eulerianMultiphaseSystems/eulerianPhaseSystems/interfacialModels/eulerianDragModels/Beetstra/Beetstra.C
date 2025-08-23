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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Beetstra.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(Beetstra, 0);
    addToRunTimeSelectionTable(eulerianDragModel, Beetstra, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::Beetstra::Beetstra
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict.lookup("residualRe"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::Beetstra::~Beetstra()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianDragModels::Beetstra::CdRe() const
{
    volScalarField alpha1
    (
        max(pair_.dispersed().volFrac(), pair_.continuous().residualAlpha())
    );

    volScalarField alpha2
    (
        max
        (
            scalar(1) - pair_.dispersed().volFrac(),
            pair_.continuous().residualAlpha()
        )
    );

    volScalarField Re(pair_.Re());
    volScalarField ReLim
    (
        "ReLim",
        max(Re, residualRe_)
    );

    volScalarField F0
    (
        "F0",
        10.0*alpha1/sqr(alpha2) + sqr(alpha2)*(1.0 + 1.5*sqrt(alpha1))
    );

    volScalarField F1
    (
        "F1",
        0.413*Re/(24.0*sqr(alpha2))*(1.0/alpha2
        + 3.0*alpha1*alpha2 + 8.4*pow(ReLim, -0.343))
        /(1.0 + pow(10.0, 3*alpha1)*pow(ReLim, -(1.0 + 4.0*alpha1)/2.0))
    );

    return 24.0*alpha2*(F0 + F1);
}


// ************************************************************************* //

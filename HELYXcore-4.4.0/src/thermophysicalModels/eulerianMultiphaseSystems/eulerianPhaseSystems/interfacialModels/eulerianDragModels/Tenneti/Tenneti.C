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

#include "Tenneti.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../SchillerNaumann/SchillerNaumann.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(Tenneti, 0);
    addToRunTimeSelectionTable(eulerianDragModel, Tenneti, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::Tenneti::Tenneti
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    SchillerNaumann_
    (
        new SchillerNaumann
        (
            dict,
            pair,
            false
        )
    ),
    residualRe_("residualRe", dimless, dict.lookup("residualRe"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::Tenneti::~Tenneti()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianDragModels::Tenneti::CdRe() const
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

    volScalarField F0
    (
        5.81*alpha1/pow3(alpha2) + 0.48*pow(alpha1, 1.0/3.0)/pow4(alpha2)
    );

    volScalarField F1
    (
        pow(alpha1, 3)*max(pair_.Re(), residualRe_)
       *(0.95 + 0.61*pow3(alpha1)/sqr(alpha2))
    );

    // Tenneti et al. correlation includes the mean pressure drag.
    // This was removed here by multiplying F by alpha2 for consistency with
    // the formulation used in OpenFOAM
    return
        SchillerNaumann_->CdRe()/(alpha2*max(pair_.Re(), residualRe_)) +
        24.0*sqr(alpha2)*(F0 + F1);
}


// ************************************************************************* //

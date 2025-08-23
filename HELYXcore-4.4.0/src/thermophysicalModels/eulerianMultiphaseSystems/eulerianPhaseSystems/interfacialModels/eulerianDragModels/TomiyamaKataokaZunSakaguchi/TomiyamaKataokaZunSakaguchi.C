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

#include "TomiyamaKataokaZunSakaguchi.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(TomiyamaKataokaZunSakaguchi, 0);
    addToRunTimeSelectionTable
    (
        eulerianDragModel,
        TomiyamaKataokaZunSakaguchi,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::TomiyamaKataokaZunSakaguchi::TomiyamaKataokaZunSakaguchi
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    residualEo_("residualEo", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::TomiyamaKataokaZunSakaguchi::~TomiyamaKataokaZunSakaguchi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianDragModels::TomiyamaKataokaZunSakaguchi::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(max(pair_.Eo(), residualEo_));

    return
        max
        (
            24.0*(1.0 + 0.15*pow(Re, 0.687))/max(Re, residualRe_),
            8.0*Eo/(3.0*(Eo + 4.0))
        )
       *max(pair_.Re(), residualRe_);
}


// ************************************************************************* //

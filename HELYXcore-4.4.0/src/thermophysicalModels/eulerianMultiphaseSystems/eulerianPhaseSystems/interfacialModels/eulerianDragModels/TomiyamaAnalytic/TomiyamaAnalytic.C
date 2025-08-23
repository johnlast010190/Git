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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "TomiyamaAnalytic.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(TomiyamaAnalytic, 0);
    addToRunTimeSelectionTable(eulerianDragModel, TomiyamaAnalytic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::TomiyamaAnalytic::TomiyamaAnalytic
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    residualEo_("residualEo", dimless, dict),
    residualE_("residualE", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::TomiyamaAnalytic::~TomiyamaAnalytic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianDragModels::TomiyamaAnalytic::CdRe() const
{
    volScalarField Eo(max(pair_.Eo(), residualEo_));
    volScalarField E(max(pair_.E(), residualE_));

    volScalarField OmEsq(max(scalar(1) - sqr(E), sqr(residualE_)));
    volScalarField rtOmEsq(sqrt(OmEsq));

    volScalarField F(max(asin(rtOmEsq) - E*rtOmEsq, residualE_)/OmEsq);

    return
        (8.0/3.0)
       *Eo
       /(
            Eo*pow(E, 2.0/3.0)/OmEsq
          + 16*pow(E, 4.0/3.0)
        )
       /sqr(F)
       *max(pair_.Re(), residualRe_);
}


// ************************************************************************* //

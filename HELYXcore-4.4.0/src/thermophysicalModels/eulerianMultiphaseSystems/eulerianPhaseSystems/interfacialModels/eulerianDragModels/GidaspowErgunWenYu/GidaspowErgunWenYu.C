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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "GidaspowErgunWenYu.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../Ergun/Ergun.H"
#include "../WenYu/WenYu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);
    addToRunTimeSelectionTable(eulerianDragModel, GidaspowErgunWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    Ergun_
    (
        new Ergun
        (
            dict,
            pair,
            false
        )
    ),
    WenYu_
    (
        new WenYu
        (
            dict,
            pair,
            false
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianDragModels::GidaspowErgunWenYu::CdRe() const
{
    return
        pos0(pair_.continuous().volFrac() - 0.8)*WenYu_->CdRe()
      + neg(pair_.continuous().volFrac() - 0.8)*Ergun_->CdRe();
}


// ************************************************************************* //

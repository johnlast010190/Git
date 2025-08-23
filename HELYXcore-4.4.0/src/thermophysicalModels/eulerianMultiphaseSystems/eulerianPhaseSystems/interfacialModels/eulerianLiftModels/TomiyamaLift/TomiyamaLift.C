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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "TomiyamaLift.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianLiftModels
{
    defineTypeNameAndDebug(TomiyamaLift, 0);
    addToRunTimeSelectionTable(eulerianLiftModel, TomiyamaLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::TomiyamaLift::TomiyamaLift
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianLiftModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::TomiyamaLift::~TomiyamaLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianLiftModels::TomiyamaLift::Cl() const
{
    volScalarField EoH(pair_.EoH2());

    volScalarField f
    (
        0.0010422*pow3(EoH) - 0.0159*sqr(EoH) - 0.0204*EoH + 0.474
    );

    return
        neg(EoH - scalar(4))*min(0.288*tanh(0.121*pair_.Re()), f)
      + pos0(EoH - scalar(4))*neg(EoH - scalar(10.7))*f
      + pos0(EoH - scalar(10.7))*(-0.288);
}


// ************************************************************************* //

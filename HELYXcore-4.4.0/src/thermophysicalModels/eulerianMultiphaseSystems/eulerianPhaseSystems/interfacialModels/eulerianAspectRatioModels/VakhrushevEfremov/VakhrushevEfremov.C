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

#include "VakhrushevEfremov.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianAspectRatioModels
{
    defineTypeNameAndDebug(VakhrushevEfremov, 0);
    addToRunTimeSelectionTable
    (
        eulerianAspectRatioModel,
        VakhrushevEfremov,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianAspectRatioModels::VakhrushevEfremov::VakhrushevEfremov
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianAspectRatioModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianAspectRatioModels::VakhrushevEfremov::~VakhrushevEfremov()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianAspectRatioModels::VakhrushevEfremov::E() const
{
    volScalarField Ta(pair_.Ta());

    return
        neg(Ta - scalar(1))*scalar(1)
      + pos0(Ta - scalar(1))*neg(Ta - scalar(39.8))
       *pow3(0.81 + 0.206*tanh(1.6 - 2*log10(max(Ta, scalar(1)))))
      + pos0(Ta - scalar(39.8))*0.24;
}


// ************************************************************************* //

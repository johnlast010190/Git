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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "Frank.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianWallLubricationModels
{
    defineTypeNameAndDebug(Frank, 0);
    addToRunTimeSelectionTable
    (
        eulerianWallLubricationModel,
        Frank,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModels::Frank::Frank
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianWallLubricationModel(dict, pair),
    Cwd_("Cwd", dimless, dict),
    Cwc_("Cwc", dimless, dict),
    p_(readScalar(dict.lookup("p")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModels::Frank::~Frank()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::eulerianWallLubricationModels::Frank::Fi() const
{
    volVectorField Ur(pair_.Ur());

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());

    volScalarField Eo(pair_.Eo());
    volScalarField yTilde(y/(Cwc_*pair_.dispersed().d()));

    return zeroGradWalls
    (
        (
            pos0(Eo - 1.0)*neg(Eo - 5.0)*exp(-0.933*Eo + 0.179)
          + pos0(Eo - 5.0)*neg(Eo - 33.0)*(0.00599*Eo - 0.0187)
          + pos0(Eo - 33.0)*0.179
        )
       *max
        (
            dimensionedScalar(dimless/dimLength, 0),
            (1.0 - yTilde)/(Cwd_*y*pow(yTilde, p_ - 1.0))
        )
       *pair_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //

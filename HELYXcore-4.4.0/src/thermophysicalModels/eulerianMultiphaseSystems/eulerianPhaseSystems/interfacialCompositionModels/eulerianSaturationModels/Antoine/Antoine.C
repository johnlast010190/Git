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

#include "Antoine.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianSaturationModels
{
    defineTypeNameAndDebug(Antoine, 0);
    addToRunTimeSelectionTable(eulerianSaturationModel, Antoine, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianSaturationModels::Antoine::Antoine(const dictionary& dict)
:
    eulerianSaturationModel(),
    A_("A", dimless, dict),
    B_("B", dimTemperature, dict),
    C_("C", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianSaturationModels::Antoine::~Antoine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::Antoine::pSat
(
    const volScalarField& T
) const
{
    return
        dimensionedScalar(dimPressure, 1)
       *exp(A_ + B_/(C_ + T));
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::Antoine::pSatPrime
(
    const volScalarField& T
) const
{
    return - pSat(T)*B_/sqr(C_ + T);
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::Antoine::lnPSat
(
    const volScalarField& T
) const
{
    return A_ + B_/(C_ + T);
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::Antoine::Tsat
(
    const volScalarField& p
) const
{
    return
        B_/(log(p*dimensionedScalar(dimless/dimPressure, 1)) - A_)
      - C_;
}


// ************************************************************************* //

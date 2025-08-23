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
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "constantSaturationConditions.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianSaturationModels
{
    defineTypeNameAndDebug(constantSaturationConditions, 0);
    addToRunTimeSelectionTable
    (
        eulerianSaturationModel,
        constantSaturationConditions,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianSaturationModels::constantSaturationConditions::
constantSaturationConditions(const dictionary& dict)
:
    eulerianSaturationModel(),
    pSat_("pSat", dimPressure, dict),
    Tsat_("Tsat", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianSaturationModels::constantSaturationConditions::
~constantSaturationConditions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::constantSaturationConditions::pSat
(
    const volScalarField& T
) const
{
    return volScalarField::New("pSat", T.mesh(), pSat_);
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::constantSaturationConditions::pSatPrime
(
    const volScalarField& T
) const
{
    return volScalarField::New
    (
        "pSatPrime",
        T.mesh(),
        dimensionedScalar(dimPressure/dimTemperature, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::constantSaturationConditions::lnPSat
(
    const volScalarField& T
) const
{
    return volScalarField::New
    (
        "lnPSat",
        T.mesh(),
        dimensionedScalar(dimless, log(pSat_.value()))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianSaturationModels::constantSaturationConditions::Tsat
(
    const volScalarField& p
) const
{
    return volScalarField::New("Tsat", p.mesh(), Tsat_);
}


// ************************************************************************* //

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
    (c) 2015-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionModels/multiBandAbsorptionEmission/multiBandAbsorptionEmission.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(multiBandAbsorptionEmission,
    0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        multiBandAbsorptionEmission,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandAbsorptionEmission::
multiBandAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeNameDerived + "Coeffs")),
    absCoeffs_(maxBands_),
    emiCoeffs_(maxBands_),
    nBands_(0)
{
    coeffsDict_.readEntry("absorptivity", absCoeffs_);
    coeffsDict_.readEntry("emissivity", emiCoeffs_);
    nBands_ = absCoeffs_.size();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandAbsorptionEmission::
~multiBandAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandAbsorptionEmission::
aCont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "aCont",
        mesh(),
        dimensionedScalar("a", dimless/dimLength, absCoeffs_[bandI])
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandAbsorptionEmission::
eCont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "eCont",
        mesh(),
        dimensionedScalar("e", dimless/dimLength, emiCoeffs_[bandI])
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandAbsorptionEmission::
ECont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "ECont",
        mesh(),
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), Zero)
    );
}


// ************************************************************************* //

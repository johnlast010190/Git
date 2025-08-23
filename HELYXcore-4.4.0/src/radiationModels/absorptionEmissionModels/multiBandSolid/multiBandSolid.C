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
    (c) 2015 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionModels/multiBandSolid/multiBandSolid.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(multiBandSolid, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        multiBandSolid,
        dictionary
    );

    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        multiBandSolid,
        dictionary,
        multiBandSolidAbsorptionEmission
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandSolid::
multiBandSolid
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
    absCoeffs_ = coeffsDict_.lookup<scalarList>("absorptivity");
    emiCoeffs_ = coeffsDict_.lookup<scalarList>("emissivity");
    nBands_ = absCoeffs_.size();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandSolid::
~multiBandSolid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandSolid::aCont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "aCont",
        mesh(),
        dimensionedScalar(dimless/dimLength, absCoeffs_[bandI])
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandSolid::eCont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "eCont",
        mesh(),
        dimensionedScalar(dimless/dimLength, emiCoeffs_[bandI])
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandSolid::ECont
(
    const label bandI
) const
{
    return absorptionEmissionModel::ECont(bandI);
}


// ************************************************************************* //

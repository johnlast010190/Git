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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/addOns/radiation/absorptionEmission/cloudAbsorptionEmission/cloudAbsorptionEmission.H"
#include "clouds/baseClasses/thermoCloud/thermoCloud.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(cloudAbsorptionEmission, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        cloudAbsorptionEmission,
        dictionary
    );
    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        cloudAbsorptionEmission,
        dictionary,
        cloudAbsorptionEmission
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::cloudAbsorptionEmission::
cloudAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& derivedClassName
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(derivedClassName + "Coeffs")),
    cloudNames_(coeffsDict_.lookup("cloudNames"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::cloudAbsorptionEmission::
~cloudAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloudAbsorptionEmission::aDisp
(
    const label
) const
{
    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "aDisp",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    forAll(cloudNames_, i)
    {
        const thermoCloud& tc
        (
            mesh_.objectRegistry::lookupObject<thermoCloud>(cloudNames_[i])
        );

        ta.ref() += tc.ap();
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloudAbsorptionEmission::eDisp
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "eDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloudAbsorptionEmission::EDisp
(
    const label bandI
) const
{
    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "EDisp",
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        )
    );

    forAll(cloudNames_, i)
    {
        const thermoCloud& tc
        (
            mesh_.objectRegistry::lookupObject<thermoCloud>(cloudNames_[i])
        );

        tE.ref() += tc.Ep();
    }

    // Total emission is 4 times the projected emission
    return 4*tE;
}


// ************************************************************************* //

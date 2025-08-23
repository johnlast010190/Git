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
    (c) 2017 Engys Ltd.
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionModels/constantAbsorptionEmission/constantAbsorptionEmission.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(constantAbsorptionEmission, 0);
    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        constantAbsorptionEmission,
        dictionary
    );
    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        constantAbsorptionEmission,
        dictionary,
        constantAbsorptionEmission
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::constantAbsorptionEmission::
constantAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeNameDerived + "Coeffs")),
    a_("absorptivity", dimless/dimLength, coeffsDict_),
    e_("emissivity", dimless/dimLength, coeffsDict_),
    E_("E", dimMass/dimLength/pow3(dimTime), coeffsDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::constantAbsorptionEmission::
~constantAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::constantAbsorptionEmission::
aCont(const label bandI) const
{
    return volScalarField::New("aCont", mesh_, a_);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::constantAbsorptionEmission::
eCont(const label bandI) const
{
    return volScalarField::New("eCont", mesh_, e_);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::constantAbsorptionEmission::
ECont(const label bandI) const
{
    return volScalarField::New("ECont", mesh_, E_);
}


// ************************************************************************* //

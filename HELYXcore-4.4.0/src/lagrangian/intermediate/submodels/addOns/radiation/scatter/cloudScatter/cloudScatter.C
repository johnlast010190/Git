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

#include "submodels/addOns/radiation/scatter/cloudScatter/cloudScatter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "clouds/baseClasses/thermoCloud/thermoCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace scatterModels
{
    defineTypeNameAndDebug(cloudScatter, 0);

    addToRunTimeSelectionTable
    (
        scatterModel,
        cloudScatter,
        dictionary
    );
    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        scatterModel,
        cloudScatter,
        dictionary,
        cloudScatter
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::scatterModels::cloudScatter::cloudScatter
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    scatterModel(dict, mesh),
    coeffsDict_(dict.subDict(typeNameDerived + "Coeffs")),
    cloudNames_(coeffsDict_.lookup("cloudNames"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::scatterModels::cloudScatter::~cloudScatter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::scatterModels::cloudScatter::sigmaEff() const
{
    tmp<volScalarField> tsigma
    (
        volScalarField::New
        (
            "sigma",
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

        tsigma.ref() += tc.sigmap();
    }

    return 3.0*tsigma;
}


// ************************************************************************* //

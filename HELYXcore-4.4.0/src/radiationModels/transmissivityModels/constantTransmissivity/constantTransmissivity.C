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

#include "transmissivityModels/constantTransmissivity/constantTransmissivity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace transmissivityModels
{
    defineTypeNameAndDebug(constantTransmissivity, 0);

    addToRunTimeSelectionTable
    (
        transmissivityModel,
        constantTransmissivity,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::transmissivityModels::constantTransmissivity::
constantTransmissivity
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    transmissivityModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    tau_(coeffsDict_.lookup<scalar>("transmissivity"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::transmissivityModels::constantTransmissivity::
~constantTransmissivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::transmissivityModels::constantTransmissivity::tauEff
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "tau",
        mesh_,
        dimensionedScalar(dimless/dimLength, tau_)
    );
}


// ************************************************************************* //

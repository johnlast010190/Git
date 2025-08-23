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

#include "radiationModels/absorptionEmissionModels/wideBandCombustion/wideBandCombustion.H"
#include "combustionModel/combustionModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(wideBandCombustion, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        wideBandCombustion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::
wideBandCombustion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    wideBand(dict, mesh, typeNameDerived)
{
    label bandi = 0;
    forAllConstIter(dictionary, coeffsDict_, iter)
    {
        if (!iter().isDict()) continue;

        iter().dict().lookup("EhrrCoeff") >> iEhrrCoeffs_[bandi];

        ++ bandi;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::~wideBandCombustion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::ECont
(
    const label bandi
) const
{
    tmp<volScalarField> E = wideBand::ECont(bandi);

    const word& name = combustionModel::combustionPropertiesName;
    E.ref() +=
        iEhrrCoeffs_[bandi]
       *mesh_.lookupObject<combustionModel>(name).Qdot()
       *(iBands_[bandi][1] - iBands_[bandi][0])
       /totalWaveLength_;

    return E;
}


// ************************************************************************* //

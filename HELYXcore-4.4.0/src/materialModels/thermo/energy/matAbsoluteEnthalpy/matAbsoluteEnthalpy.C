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
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "matAbsoluteEnthalpy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(matAbsoluteEnthalpy, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        matAbsoluteEnthalpy,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matAbsoluteEnthalpy::matAbsoluteEnthalpy
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
}


Foam::autoPtr<Foam::matAbsoluteEnthalpy>
Foam::matAbsoluteEnthalpy::clone() const
{
    return autoPtr<matAbsoluteEnthalpy>(new matAbsoluteEnthalpy(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::matAbsoluteEnthalpy::updateTable(const word& modelName)
{
    // Model names
    const wordList modList({CpvModel::typeName, heModel::typeName});
    dep_.setSize(modList.size());

    // Required dependencies
    const wordList depNames({CpModel::typeName, haModel::typeName});

    // Dependency indices
    const labelList depInds({Cp, ha});

    // Create up links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::matAbsoluteEnthalpy::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, Cpv)
    castMaterial(modelName, he)

    return nullptr;
}


// ************************************************************************* //

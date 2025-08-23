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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pvInvertFunc.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pvInvertFunc, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        pvInvertFunc,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pvInvertFunc::pvInvertFunc
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
    pvInvertFunc::read();
}


Foam::autoPtr<Foam::pvInvertFunc>
Foam::pvInvertFunc::clone() const
{
    return autoPtr<pvInvertFunc>(new pvInvertFunc(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pvInvertFunc::updateTable(const word& modelName)
{
    // Model names
    const wordList modList({pvInvertModel::typeName});
    dep_.setSize(modList.size());

    // Required dependencies
    const wordList depNames(modList.size(), pvModel::typeName);

    // Dependency indices
    const labelList depInds(modList.size(), pvInd);

    // Create up links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::pvInvertFunc::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, pvInvert)
    return nullptr;
}


bool Foam::pvInvertFunc::read()
{
    initialise("p");
    initialise("T");

    // Make sure that properties are updated
    materialTables_.readScalarProperties();
    Pc_ = materialTables_.characterProperties("Pc");
    Tc_ = materialTables_.characterProperties("Tc");
    Pt_ = materialTables_.characterProperties("Pt");
    Tt_ = materialTables_.characterProperties("Tt");
    Tb_ = materialTables_.characterProperties("Tb");

    return true;
}


// ************************************************************************* //

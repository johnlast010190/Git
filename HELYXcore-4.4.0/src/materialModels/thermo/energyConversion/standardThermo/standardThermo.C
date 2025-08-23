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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "standardThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(standardThermo, 0);
    addToRunTimeSelectionTable(materialModel, standardThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::standardThermo::standardThermo
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
    standardThermo::read();
}


Foam::autoPtr<Foam::standardThermo> Foam::standardThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<standardThermo>
    (
        new standardThermo(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standardThermo::updateTable(const word& modelName)
{
    // Model names
    const wordList modList
    ({
        gammaModel::typeName,
        gModel::typeName,
        aModel::typeName
    });
    dep_.setSize(modList.size());

    // Required dependencies
    const List<wordList> depNames
    ({
        {CpModel::typeName, CpMCvModel::typeName},
        {haModel::typeName, sModel::typeName},
        {eaModel::typeName, sModel::typeName}
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {Cp, CpMCv},
        {ha, s},
        {eaInd, s}
    });

    // Create up links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::standardThermo::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, gamma)
    castMaterial(modelName, g)
    castMaterial(modelName, a)

    return nullptr;
}


bool Foam::standardThermo::read()
{
    initialise("p");
    initialise("T");
    return true;
}


// ************************************************************************* //

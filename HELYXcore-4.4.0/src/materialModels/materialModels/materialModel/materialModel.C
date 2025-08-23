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
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "materialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialModel, 0);
    defineRunTimeSelectionTable(materialModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModel::materialModel
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    regIOobject
    (
        IOobject
        (
            name,
            obr.time().timeName(),
            obr.subRegistry("materialModels", true, false)
        )
    ),
    matObr_(obr.subRegistry("materialModels")),
    materialTables_(matObr_.lookupObjectRef<materialTables>("materialTables")),
    phaseName_(phaseName),
    specieName_(specieName),
    name_(name),
    mesh_(fvSolutionRegistry::getMesh(obr)),
    obr_(obr),
    dict_(&dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModel::~materialModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::materialModel::updateDictPtr()
{
    const bool isPhaseMixture = (phaseName_ == word::null) ? false : true;
    const bool isSpecieMixture = (specieName_ == word::null) ? false : true;

    const dictionary& matDict = materialsDict();

    const auto i = name_.rfind('.');

    word coeffDictName;

    if (i == std::string::npos || i == 0)
    {
        coeffDictName = name_;
    }
    else
    {
        coeffDictName = name_.substr(0, i);
    }
    coeffDictName += "Coeffs";

    // Model dictionary
    dict_ =
        isPhaseMixture ?
        (
            isSpecieMixture
          ? &matDict.subDict(phaseName_).subDict(specieName_).optionalSubDict(coeffDictName)
          : &matDict.subDict(phaseName_).optionalSubDict(coeffDictName)
        )
      : (
            isSpecieMixture
          ? &matDict.subDict(specieName_).optionalSubDict(coeffDictName)
          : &matDict.optionalSubDict(coeffDictName)
        );

    if (dict_ == nullptr)
    {
        FatalErrorInFunction
            << "Sub-dictionary not found."
            << exit(FatalError);
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const word& modelName,
    const wordList& requiredModels,
    const labelList& enumInds,
    const label depIndex
)
{
    const boolList compulsory(enumInds.size(), true);
    fill(checkModel, modelName, requiredModels, enumInds, compulsory, depIndex);
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const word& modelName,
    const word& requiredModel,
    const label enumInd,
    const label depNumber
)
{
    const labelList enumInds({enumInd});
    const wordList requiredModels({requiredModel});
    fill(checkModel, modelName, requiredModels, enumInds, depNumber);
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const word& modelName,
    const wordList& modelNames,
    const labelList& enumInds,
    const boolList& compulsory,
    const label depNumber
)
{
    if (checkModel == modelName)
    {
        boolList modelsAdeed(modelNames.size(), false);
        const word tName(materialTables_.tableName(phaseName_, specieName_));

        depSubList oneDep;
        label size = 0;
        if (materialTables_.sTable().found(tName))
        {
            const matScalarTable& sModels =
                materialTables_.sTable(phaseName_, specieName_);
            if (sModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (sModels.found(name))
                    {
                        oneDep.model = sModels[modelName];
                        oneDep.dependencies.setSize(++size);
                        modelsAdeed[i] = true;
                        sMod_.set(enumInds[i], sModels[name]);
                        oneDep.dependencies.set(size - 1, sModels[name]);
                    }
                }
            }
        }

        if (materialTables_.vTable().found(tName))
        {
            const matVectorTable& vModels =
                materialTables_.vTable(phaseName_, specieName_);
            if (vModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (vModels.found(name))
                    {
                        oneDep.model = vModels[modelName];
                        oneDep.dependencies.setSize(++size);
                        modelsAdeed[i] = true;
                        vMod_.set(enumInds[i], vModels[name]);
                        oneDep.dependencies.set(size - 1, vModels[name]);
                    }
                }
            }
        }

        if (materialTables_.tTable().found(tName))
        {
            const matTensorTable& tModels =
                materialTables_.tTable(phaseName_, specieName_);
            if (tModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (tModels.found(name))
                    {
                        oneDep.model = tModels[modelName];
                        oneDep.dependencies.setSize(++size);
                        modelsAdeed[i] = true;
                        tMod_.set(enumInds[i], tModels[name]);
                        oneDep.dependencies.set(size - 1, tModels[name]);
                    }
                }
            }
        }
        if (oneDep.model != nullptr)
        {
            dep_[depNumber] = oneDep;
        }
        bool allUpdated = true;
        forAll(modelsAdeed, i)
        {
            allUpdated =
                (compulsory[i] == true && modelsAdeed[i] == false)
              ? false
              : allUpdated;
        }
        if (!allUpdated)
        {
            FatalErrorInFunction
                << "At least one of the required models for \""
                << modelName << "\""
                << " from the list: " << nl << modelNames << " not found."
                << exit(FatalError);
        }
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const wordList& modList,
    const wordList& dependencies,
    const labelList& enumInds
)
{
    forAll(modList, i)
    {
        fill(checkModel, modList[i], dependencies[i], enumInds[i], i);
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const wordList& modList,
    const List<wordList>& dependencies,
    const List<labelList>& enumInds
)
{
    forAll(modList, i)
    {
        fill(checkModel, modList[i], dependencies[i], enumInds[i], i);
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const wordList& modList,
    const wordList& dependencies,
    const labelList& enumInds,
    const boolList& compulsory
)
{
    forAll(modList, i)
    {
        fill
        (
            checkModel,
            modList[i],
            {dependencies[i]},
            {enumInds[i]},
            {compulsory[i]},
            i
        );
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const wordList& modList,
    const List<wordList>& dependencies,
    const List<labelList>& enumInds,
    const List<boolList>& compulsory
)
{
    forAll(modList, i)
    {
        fill
        (
            checkModel,
            modList[i],
            dependencies[i],
            enumInds[i],
            compulsory[i],
            i
        );
    }
}


const Foam::dictionary& Foam::materialModel::materialsDict() const
{
    return materialTables_.dict();
}


// ************************************************************************* //

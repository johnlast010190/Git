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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "materialModels/materialTables/materialTables.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coefficientMixture> Foam::coefficientMixture::New
(
    const objectRegistry& obr,
    const scalarField& stoichCoeffs,
    const wordList& speciesNames,
    const word& phaseName,
    const word& name
)
{
    const objectRegistry& matObr = obr.subRegistry("materialModels");
    materialTables& mat =
        matObr.lookupObjectRef<materialTables>("materialTables");
    bool consistentModels = true;

    // Coefficient mixing part
    UPtrList<baseModels<scalar>> modelsForCoeffMix(stoichCoeffs.size());
    scalarField chemistryCoffs(stoichCoeffs.size());
    scalarField molarWeights(stoichCoeffs.size());

    // Initialise name of the model for the first specie
    matScalarTable& speciesTable = mat.sTable(phaseName, speciesNames[0]);
    const word haModel(speciesTable[haModel::typeName]->type());
    const word sModel(speciesTable[sModel::typeName]->type());

    forAll(modelsForCoeffMix, modeli)
    {
        const word specieName = speciesNames[modeli];
        matScalarTable& speciesTable = mat.sTable(phaseName, specieName);
        const scalar W = (*speciesTable[WModel::typeName])[0];

        // Check model types
        const word haModeli(speciesTable[haModel::typeName]->type());
        const word sModeli(speciesTable[sModel::typeName]->type());

        if (haModel == haModeli && sModel == sModeli)
        {
            modelsForCoeffMix.set(modeli, speciesTable[haModel::typeName]);
            chemistryCoffs[modeli] = W*stoichCoeffs[modeli];
        }
        else
        {
            consistentModels = false;
            break;
        }
    }

    // Consturct only models that have have coefficient mixtures and
    // if they have all the models of the same type in the chemical equation.
    typename dictionaryConstructorTable::iterator iter =
        dictionaryConstructorTable_().find(haModel);

    if (consistentModels && iter != dictionaryConstructorTable_().end())
    {
        const auto ctor =
            ctorTableLookup
            (
                "material type",
                dictionaryConstructorTable_(),
                haModel
            );
        return
            autoPtr<coefficientMixture>
            (
                ctor(modelsForCoeffMix, chemistryCoffs, name)
            );
    }
    else
    {
        return autoPtr<coefficientMixture>();
    }
}


// ************************************************************************* //

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

#include "homogeneousThermoMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(homogeneousThermoMixture, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        homogeneousThermoMixture,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::speciesMassFractions&
Foam::homogeneousThermoMixture::lookupOrConstructBase()
{
    const word name("massVolumeFractions");
    if (!obr_.foundObject<speciesMassFractions>(name))
    {
        obr_.store
        (
            new speciesMassFractions
            (
                materialsDict(),
                wordList({"b"}),
                obr_,
                word::null
            )
        );
    }
    return obr_.lookupObject<speciesMassFractions>(name);
}


Foam::homogeneousThermoMixture::homogeneousThermoMixture
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    frac_(lookupOrConstructBase()),
    b_(frac_.Y("b"))
{}


Foam::autoPtr<Foam::homogeneousThermoMixture>
Foam::homogeneousThermoMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<homogeneousThermoMixture>
    (
        new homogeneousThermoMixture(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::homogeneousThermoMixture::updateTable(const word& modelName)
{
    const wordList speciesNames = dict_->lookup<wordList>("species");
    const word ttName = materialTables_.tableName(phaseName_, word::null);

    const HashTable<matScalarTable>& sModels = materialTables_.sTable();
    word funcName = IOobject::group();

    if (IOobject::name().find("Model") != string::npos)
    {
        funcName = string(this->name()).replace("Model", "");
    }
    const auto i = funcName.find('.');
    if (i != 0)
    {
        funcName = funcName.substr(0, i);
    }
    if (materialTables_.foundModel<scalar>(ttName, funcName))
    {
        sMod_.setSize(speciesNames.size());
        dep_.setSize(1);
        dep_[0].model = sModels[ttName][funcName];
        dep_[0].dependencies.setSize(speciesNames.size());

        word tName =
            materialTables_.tableName(phaseName_, "reactants");
        sMod_.set(reactantsInd, sModels[tName][funcName]);
        dep_[0].dependencies.set(reactantsInd, sModels[tName][funcName]);

        tName = materialTables_.tableName(phaseName_, "products");
        sMod_.set(productsInd, sModels[tName][funcName]);
        dep_[0].dependencies.set(productsInd, sModels[tName][funcName]);
        harmonic_ =
            (
                funcName == rhoModel::typeName
             || funcName == psiModel::typeName
             || funcName == WModel::typeName
             || funcName == "buoyancy"
            );
    }
}


Foam::baseModels<Foam::scalar>* Foam::homogeneousThermoMixture::castScalarModel
(
    const word& modelName
)
{
    const word mixType =
        dict_->lookupOrDefault<word>("mixture", "homogeneousThermoMixture");

    if (homogeneousReactantMixtureModel::typeName == mixType)
    {
        return dynamic_cast<homogeneousReactantMixtureModel*>(this);
    }
    else if (homogeneousProductMixtureModel::typeName == mixType)
    {
        return dynamic_cast<homogeneousProductMixtureModel*>(this);
    }
    else if (homogeneousThermoMixtureModel::typeName == mixType)
    {
        return dynamic_cast<homogeneousThermoMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::homogeneousThermoMixture::homogeneousThermoMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[reactantsInd].funcType(),
            obr_,
            frac_.fractions()[0].mesh(),
            dimensionedScalar(sMod_[reactantsInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = homogeneousThermoMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            homogeneousThermoMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::homogeneousThermoMixture::homogeneousThermoMixtureInternal() const
{
    const scalarField reactants(sMod_[reactantsInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& b = b_.primitiveField();
    tmp<scalarField> tMix(new scalarField(reactants.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture(b[celli], reactants[celli], products[celli]);
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::homogeneousThermoMixture::homogeneousThermoMixturePatch
(
    const label patchi
) const
{
    const scalarField reactants(sMod_[reactantsInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& b = b_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(reactants.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                b[patchFacei],
                reactants[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar Foam::homogeneousThermoMixture::homogeneousThermoMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        b_[celli],
        sMod_[reactantsInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar Foam::homogeneousThermoMixture::homogeneousThermoMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


bool Foam::homogeneousThermoMixture::read()
{
    return true;
}


// ************************************************************************* //

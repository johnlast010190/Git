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

#include "egrThermoMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(egrThermoMixture, 0);
    addToRunTimeSelectionTable(materialModel, egrThermoMixture, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::speciesMassFractions&
Foam::egrThermoMixture::lookupOrConstructBase()
{
    const word name("massVolumeFractions");
    if (!obr_.foundObject<speciesMassFractions>(name))
    {
        obr_.store
        (
            new speciesMassFractions
            (
                materialsDict(),
                wordList({"ft", "b", "egr"}),
                obr_,
                word::null
            )
        );
    }
    return obr_.lookupObject<speciesMassFractions>(name);
}


Foam::egrThermoMixture::egrThermoMixture
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
    ft_(frac_.Y("ft")),
    b_(frac_.Y("b")),
    egr_(frac_.Y("egr")),
    stoicRatio_
    (
        materialsDict().lookup<scalar>("stoichiometricAirFuelMassRatio")
    )
{}


Foam::autoPtr<Foam::egrThermoMixture> Foam::egrThermoMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<egrThermoMixture>
    (
        new egrThermoMixture(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::egrThermoMixture::updateTable(const word& modelName)
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
            materialTables_.tableName(phaseName_, "fuel");
        sMod_.set(fuelInd, sModels[tName][funcName]);
        dep_[0].dependencies.set(fuelInd, sModels[tName][funcName]);

        tName =
            materialTables_.tableName(phaseName_, "oxidant");
        sMod_.set(oxidantInd, sModels[tName][funcName]);
        dep_[0].dependencies.set(oxidantInd, sModels[tName][funcName]);

        tName =
            materialTables_.tableName(phaseName_, "burntProducts");
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


Foam::baseModels<Foam::scalar>* Foam::egrThermoMixture::castScalarModel
(
    const word& modelName
)
{
    const word mixType =
        dict_->lookupOrDefault<word>("mixture", "egrThermoMixture");

    if (mixType == egrReactantMixtureModel::typeName)
    {
        return dynamic_cast<egrReactantMixtureModel*>(this);
    }
    else if (egrProductMixtureModel::typeName == mixType)
    {
        return dynamic_cast<egrProductMixtureModel*>(this);
    }
    else if (egrThermoMixtureModel::typeName == mixType)
    {
        return dynamic_cast<egrThermoMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::egrThermoMixture::egrReactantMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[oxidantInd].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = egrReactantMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] = egrReactantMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrReactantMixtureInternal() const
{
    const scalarField oxidant(sMod_[oxidantInd].primitiveField());
    const scalarField fuel(sMod_[fuelInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& ft = ft_.primitiveField();
    const scalarField& egr = egr_.primitiveField();
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture
            (
                ft[celli],
                1,
                egr[celli],
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrReactantMixturePatch(const label patchi) const
{
    const scalarField oxidant(sMod_[oxidantInd].boundaryField()[patchi]);
    const scalarField fuel(sMod_[fuelInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& ft = ft_.boundaryField()[patchi];
    const scalarField& egr = egr_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                ft[patchFacei],
                1,
                egr[patchFacei],
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar Foam::egrThermoMixture::egrReactantMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        1,
        egr_[celli],
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar Foam::egrThermoMixture::egrReactantMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrProductMixtureInternal() const
{
    const scalarField oxidant(sMod_[oxidantInd].primitiveField());
    const scalarField fuel(sMod_[fuelInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& ft = ft_.primitiveField();
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture
            (
                ft[celli],
                0,
                0,
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrProductMixturePatch(const label patchi) const
{
    const scalarField oxidant(sMod_[oxidantInd].boundaryField()[patchi]);
    const scalarField fuel(sMod_[fuelInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& ft = ft_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                ft[patchFacei],
                0,
                0,
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::tmp<Foam::volScalarField>
Foam::egrThermoMixture::egrProductMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[oxidantInd].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = egrProductMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] = egrProductMixturePatch(patchi);
    }
    return tMixture;
}


Foam::scalar
Foam::egrThermoMixture::egrProductMixtureCell(const label celli) const
{
    return combustionMixture
    (
        ft_[celli],
        0,
        0,
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar
Foam::egrThermoMixture::egrProductMixtureValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::egrThermoMixture::egrThermoMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[oxidantInd].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = egrThermoMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] = egrThermoMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrThermoMixtureInternal() const
{
    const scalarField oxidant(sMod_[oxidantInd].primitiveField());
    const scalarField fuel(sMod_[fuelInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& ft = ft_.primitiveField();
    const scalarField& b = b_.primitiveField();
    const scalarField& egr = egr_.primitiveField();
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture
            (
                ft[celli],
                b[celli],
                egr[celli],
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::egrThermoMixture::egrThermoMixturePatch(const label patchi) const
{
    const scalarField oxidant(sMod_[oxidantInd].boundaryField()[patchi]);
    const scalarField fuel(sMod_[fuelInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& ft = ft_.boundaryField()[patchi];
    const scalarField& b = b_.boundaryField()[patchi];
    const scalarField& egr = egr_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                ft[patchFacei],
                b[patchFacei],
                egr[patchFacei],
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar
Foam::egrThermoMixture::egrThermoMixtureCell(const label celli) const
{
    return combustionMixture
    (
        ft_[celli],
        b_[celli],
        egr_[celli],
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar
Foam::egrThermoMixture::egrThermoMixtureValue(scalar p, scalar T) const
{
    NotImplemented;
}


bool Foam::egrThermoMixture::read()
{
    return true;
}


// ************************************************************************* //

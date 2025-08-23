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

#include "inhomogeneousThermoMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inhomogeneousThermoMixture, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        inhomogeneousThermoMixture,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::speciesMassFractions&
Foam::inhomogeneousThermoMixture::lookupOrConstructBase()
{
    const word name("massVolumeFractions");
    if (!obr_.foundObject<speciesMassFractions>(name))
    {
        obr_.store
        (
            new speciesMassFractions
            (
                materialsDict(),
                wordList({"ft", "b"}),
                obr_,
                word::null
            )
        );
    }
    return obr_.lookupObject<speciesMassFractions>(name);
}


Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixture
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
    stoicRatio_
    (
        materialsDict().lookup<scalar>("stoichiometricAirFuelMassRatio")
    )
{}


Foam::autoPtr<Foam::inhomogeneousThermoMixture>
Foam::inhomogeneousThermoMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<inhomogeneousThermoMixture>
    (
        new inhomogeneousThermoMixture(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inhomogeneousThermoMixture::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>*
Foam::inhomogeneousThermoMixture::castScalarModel
(
    const word& modelName
)
{
    const word mixType =
        dict_->lookupOrDefault<word>("mixture", "inhomogeneousThermoMixture");

    if (mixType == inhomogeneousReactantMixtureModel::typeName)
    {
        return dynamic_cast<inhomogeneousReactantMixtureModel*>(this);
    }
    else if (inhomogeneousProductMixtureModel::typeName == mixType)
    {
        return dynamic_cast<inhomogeneousProductMixtureModel*>(this);
    }
    else if (inhomogeneousThermoMixtureModel::typeName == mixType)
    {
        return dynamic_cast<inhomogeneousThermoMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousReactantMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[0].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = inhomogeneousReactantMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            inhomogeneousReactantMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousReactantMixtureInternal() const
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
                1,
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousReactantMixturePatch
(
    const label patchi
) const
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
                1,
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar Foam::inhomogeneousThermoMixture::inhomogeneousReactantMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        1,
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar Foam::inhomogeneousThermoMixture::inhomogeneousReactantMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousProductMixtureInternal() const
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
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousProductMixturePatch
(
    const label patchi
) const
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
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::tmp<Foam::volScalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousProductMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[0].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = inhomogeneousProductMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            inhomogeneousProductMixturePatch(patchi);
    }
    return tMixture;
}


Foam::scalar
Foam::inhomogeneousThermoMixture::inhomogeneousProductMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        0,
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar Foam::inhomogeneousThermoMixture::inhomogeneousProductMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixtureGeometric() const
{
    tmp<volScalarField> tMixture
    (
        volScalarField::New
        (
            sMod_[oxidantInd].funcType(),
            obr_,
            frac_.fractions()[0].mesh(),
            dimensionedScalar(sMod_[oxidantInd].dimensions(), 0)
        )
    );
    volScalarField& mixture = tMixture.ref();
    mixture.primitiveFieldRef() = inhomogeneousThermoMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            inhomogeneousThermoMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixtureInternal() const
{
    const scalarField oxidant(sMod_[oxidantInd].primitiveField());
    const scalarField fuel(sMod_[fuelInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& ft = ft_.primitiveField();
    const scalarField& b = b_.primitiveField();
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture
            (
                ft[celli],
                b[celli],
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixturePatch
(
    const label patchi
) const
{
    const scalarField oxidant(sMod_[oxidantInd].boundaryField()[patchi]);
    const scalarField fuel(sMod_[fuelInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& ft = ft_.boundaryField()[patchi];
    const scalarField& b = b_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                ft[patchFacei],
                b[patchFacei],
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar
Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        b_[celli],
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar
Foam::inhomogeneousThermoMixture::inhomogeneousThermoMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


bool Foam::inhomogeneousThermoMixture::read()
{
    return true;
}


// ************************************************************************* //

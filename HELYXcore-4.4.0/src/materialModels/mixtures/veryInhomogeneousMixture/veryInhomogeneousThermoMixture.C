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

#include "veryInhomogeneousThermoMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(veryInhomogeneousThermoMixture, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        veryInhomogeneousThermoMixture,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::speciesMassFractions&
Foam::veryInhomogeneousThermoMixture::lookupOrConstructBase()
{
    const word name("massVolumeFractions");
    if (!obr_.foundObject<speciesMassFractions>(name))
    {
        obr_.store
        (
            new speciesMassFractions
            (
                materialsDict(),
                wordList({"ft", "b", "veryInhomogeneous"}),
                obr_,
                word::null
            )
        );
    }
    return obr_.lookupObject<speciesMassFractions>(name);
}


Foam::veryInhomogeneousThermoMixture::veryInhomogeneousThermoMixture
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
    fu_(frac_.Y("fu")),
    b_(frac_.Y("b")),
    stoicRatio_
    (
        materialsDict().lookup<scalar>("stoichiometricAirFuelMassRatio")
    )
{}


Foam::autoPtr<Foam::veryInhomogeneousThermoMixture>
Foam::veryInhomogeneousThermoMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<veryInhomogeneousThermoMixture>
    (
        new veryInhomogeneousThermoMixture(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::veryInhomogeneousThermoMixture::updateTable(const word& modelName)
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
Foam::veryInhomogeneousThermoMixture::castScalarModel
(
    const word& modelName
)
{
    const word mixType =
        dict_->lookupOrDefault<word>
        (
            "mixture",
            "veryInhomogeneousThermoMixture"
        );

    if (mixType == veryInhomogeneousReactantMixtureModel::typeName)
    {
        return dynamic_cast<veryInhomogeneousReactantMixtureModel*>(this);
    }
    else if (veryInhomogeneousProductMixtureModel::typeName == mixType)
    {
        return dynamic_cast<veryInhomogeneousProductMixtureModel*>(this);
    }
    else if (veryInhomogeneousThermoMixtureModel::typeName == mixType)
    {
        return dynamic_cast<veryInhomogeneousThermoMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousReactantMixtureGeometric() const
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
    mixture.primitiveFieldRef() = veryInhomogeneousReactantMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            veryInhomogeneousReactantMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousReactantMixtureInternal() const
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
                ft[celli],
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousReactantMixturePatch
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
                ft[patchFacei],
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousReactantMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        ft_[celli],
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousReactantMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousProductMixtureInternal() const
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
                fres(ft[celli]),
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousProductMixturePatch
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
                fres(ft[patchFacei]),
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::tmp<Foam::volScalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousProductMixtureGeometric() const
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
    mixture.primitiveFieldRef() = veryInhomogeneousProductMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            veryInhomogeneousProductMixturePatch(patchi);
    }
    return tMixture;
}


Foam::scalar
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousProductMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        fres(ft_[celli]),
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousProductMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousThermoMixtureGeometric() const
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
    mixture.primitiveFieldRef() = veryInhomogeneousThermoMixtureInternal();

    forAll(mixture.boundaryField(), patchi)
    {
        mixture.boundaryFieldRef()[patchi] =
            veryInhomogeneousThermoMixturePatch(patchi);
    }
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::
veryInhomogeneousThermoMixtureInternal() const
{
    const scalarField oxidant(sMod_[oxidantInd].primitiveField());
    const scalarField fuel(sMod_[fuelInd].primitiveField());
    const scalarField products(sMod_[productsInd].primitiveField());
    const scalarField& ft = ft_.primitiveField();
    const scalarField& fu = fu_.primitiveField();
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, celli)
    {
        mix[celli] =
            combustionMixture
            (
                ft[celli],
                fu[celli],
                oxidant[celli],
                fuel[celli],
                products[celli]
            );
    }
    return tMix;
}


Foam::tmp<Foam::scalarField>
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousThermoMixturePatch
(
    const label patchi
) const
{
    const scalarField oxidant(sMod_[oxidantInd].boundaryField()[patchi]);
    const scalarField fuel(sMod_[fuelInd].boundaryField()[patchi]);
    const scalarField products(sMod_[productsInd].boundaryField()[patchi]);
    const scalarField& ft = ft_.boundaryField()[patchi];
    const scalarField& fu = fu_.boundaryField()[patchi];
    tmp<scalarField> tMix(new scalarField(oxidant.size(), 0));
    scalarField& mix = tMix.ref();
    forAll(mix, patchFacei)
    {
        mix[patchFacei] =
            combustionMixture
            (
                ft[patchFacei],
                fu[patchFacei],
                oxidant[patchFacei],
                fuel[patchFacei],
                products[patchFacei]
            );
    }
    return tMix;
}


Foam::scalar
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousThermoMixtureCell
(
    const label celli
) const
{
    return combustionMixture
    (
        ft_[celli],
        fu_[celli],
        sMod_[oxidantInd][celli],
        sMod_[fuelInd][celli],
        sMod_[productsInd][celli]
    );
}


Foam::scalar
Foam::veryInhomogeneousThermoMixture::veryInhomogeneousThermoMixtureValue
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


bool Foam::veryInhomogeneousThermoMixture::read()
{
    return true;
}


// ************************************************************************* //

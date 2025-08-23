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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "functionTAnIso.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionTAnIso, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        functionTAnIso,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionTAnIso::functionTAnIso
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
    vMod_.setSize(1);
    dep_.setSize(1);
    functionTAnIso::read();
}


Foam::autoPtr<Foam::functionTAnIso>
Foam::functionTAnIso::clone() const
{
    return autoPtr<functionTAnIso>(new functionTAnIso(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionTAnIso::updateTable(const word& modelName)
{
    // Note mixed fill functions are not supported yet
    // Hance this has to be done manually
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);
    const matVectorTable& vModels =
        materialTables_.vTable(phaseName_, specieName_);

    // Model names
    const word& kappaName = kappaModel::typeName;
    const word& vKappaName = vKappaModel::typeName;

    if (modelName == kappaName)
    {
        vMod_.set(0, vModels[vKappaName]);
        dep_[0].model = models[kappaName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, vModels[vKappaName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::functionTAnIso::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, kappa)

    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::functionTAnIso::castVectorModel
(
    const word& modelName
)
{
    castMaterial(modelName, vKappa)

    return nullptr;
}


Foam::scalar Foam::functionTAnIso::kappaCell(const label celli) const
{
    return mag(vMod_[0][celli]);
}


Foam::scalar Foam::functionTAnIso::kappaValue(scalar p, scalar T) const
{
    return mag(vMod_[0].value(p, T));
}


Foam::tmp<Foam::scalarField> Foam::functionTAnIso::kappaPatch
(
    const label patchi
) const
{
    return mag(vMod_[0].boundaryField()[patchi]);
}


Foam::tmp<Foam::scalarField> Foam::functionTAnIso::kappaInternal() const
{
    return mag(vMod_[0].primitiveField());
}


Foam::tmp<Foam::volScalarField>
Foam::functionTAnIso::kappaGeometric() const
{
    tmp<volScalarField> tkappa
    (
        volScalarField::New
        (
            "kappa",
            obr_,
            mesh(),
            dimensionedScalar(kappaModel::modelDims, 0)
        )
    );
    volScalarField& kap = tkappa.ref();
    kap.primitiveFieldRef() = kappaInternal();

    forAll(kap.boundaryField(), patchi)
    {
        kap.boundaryFieldRef()[patchi].forceAssign(kappaPatch(patchi));
    }

    return tkappa;
}


bool Foam::functionTAnIso::read()
{
    initialise("T");
    kappa_.reset(Function1<vector>::New("kappa", *dict_));

    return true;
}


// ************************************************************************* //

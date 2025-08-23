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
    (c) 2022-2025 Engys Ltd.
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "SolidThermo.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::SolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    matScalarTable& sModels =
        this->materials_.sTable(this->phaseName_, word::null);
    const baseModels<scalar>& TMod = (*sModels[TModel::typeName]);
    const baseModels<scalar>& rhoMod = (*sModels[rhoModel::typeName]);
    const baseModels<scalar>& kappaMod = (*sModels[kappaModel::typeName]);
    const baseModels<scalar>& heMod = (*sModels[heModel::typeName]);
    const baseModels<scalar>& CpMod = (*sModels[CpModel::typeName]);
    const baseModels<scalar>& CvMod = (*sModels[CvModel::typeName]);

    this->T_.primitiveFieldRef() = TMod.primitiveField();
    this->Cp_.primitiveFieldRef() = CpMod.primitiveField();
    this->Cv_.primitiveFieldRef() = CvMod.primitiveField();
    this->rho_.primitiveFieldRef() = rhoMod.primitiveField();
    this->kappa_.primitiveFieldRef() = kappaMod.primitiveField();

    forAll(this->T_.boundaryField(), patchi)
    {
        if (this->T_.boundaryFieldRef()[patchi].fixesValue())
        {
            this->he().boundaryFieldRef()[patchi].forceAssign
            (
                heMod.boundaryField()[patchi]
            );
        }
        else
        {
            this->T_.boundaryFieldRef()[patchi].forceAssign
            (
                TMod.boundaryField()[patchi]
            );
        }
        this->rho_.boundaryFieldRef()[patchi].forceAssign
        (
            rhoMod.boundaryField()[patchi]
        );
        this->kappa_.boundaryFieldRef()[patchi].forceAssign
        (
            kappaMod.boundaryField()[patchi]
        );
        this->Cp_.boundaryFieldRef()[patchi].forceAssign
        (
            CpMod.boundaryField()[patchi]
        );
        this->Cv_.boundaryFieldRef()[patchi].forceAssign
        (
            CvMod.boundaryField()[patchi]
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::SolidThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo<BasicSolidThermo, MixtureType>(obr, phaseName),
    coorFramePtr_(nullptr)
{
    calculate();
}


template<class BasicSolidThermo, class MixtureType>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::SolidThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo<BasicSolidThermo, MixtureType>(obr, dict, phaseName),
    coorFramePtr_(nullptr)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::~SolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
const Foam::coordinateSystem&
Foam::SolidThermo<BasicSolidThermo, MixtureType>::loadCoor() const
{
    if (!coorFramePtr_)
    {
        Info<< "Adding reference frame for anisotropic solid\n" << endl;
        const word frameName =
            this->properties().isDict("kappaModelCoeffs")
          ? this->properties().subDict("kappaModelCoeffs").template
            lookupOrDefault<word>("referenceFrame", word::null)
          : word::null;
        if (frameName != word::null)
        {
            coorFramePtr_ = &coordinateFrame::New(this->T_.mesh(), frameName);
        }
        else
        {
            Info<< "Using global reference frame for anisotropic solid\n"
                << endl;
            coorFramePtr_ = coordinateFrame::globalFrame(this->T_.mesh());
        }
    }
    return coorFramePtr_->coorSys();
}


template<class BasicSolidThermo, class MixtureType>
void Foam::SolidThermo<BasicSolidThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volVectorField>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::Kappa() const
{
    return
        this->materials_.vTable
        (
            this->phaseName_,
            word::null
        )[vKappaModel::typeName]->operator()();
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::vectorField>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::Kappa
(
    const label patchi
) const
{
    return
        this->materials_.vTable
        (
            this->phaseName_,
            word::null
        )[vKappaModel::typeName]->boundaryField()[patchi];
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volSymmTensorField>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::KappaLocal() const
{
    const fvMesh& mesh = this->T_.mesh();

    const tmp<volVectorField> tKappa(Kappa());
    const volVectorField& Kappa = tKappa();

    tmp<volSymmTensorField> tKappaLocal
    (
        volSymmTensorField::New
        (
            "KappaLocal",
            mesh,
            dimensionedSymmTensor(Kappa.dimensions(), Zero)
        )
    );
    volSymmTensorField& KappaLocal = tKappaLocal.ref();

    KappaLocal.primitiveFieldRef() =
        loadCoor().transformPrincipal(mesh.C(), Kappa);
    forAll(KappaLocal.boundaryField(), patchi)
    {
        fvPatchField<symmTensor>& pf = KappaLocal.boundaryFieldRef()[patchi];
        if (!isA<emptyFvPatch>(pf.patch()))
        {
            pf.forceAssign
            (
                loadCoor().transformPrincipal
                (
                    mesh.C().boundaryField()[patchi],
                    Kappa.boundaryField()[patchi]
                )
            );
        }
    }

    return tKappaLocal;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::symmTensorField>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::KappaLocal
(
    const label patchi
) const
{
    return loadCoor().transformPrincipal
    (
        this->T_.mesh().C().boundaryField()[patchi],
        Kappa(patchi)
    );
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volSymmTensorField>
Foam::SolidThermo<BasicSolidThermo, MixtureType>::alphahLocal() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("Anialpha", this->phaseName_),
        KappaLocal()/this->Cp()
    );
}


// ************************************************************************* //

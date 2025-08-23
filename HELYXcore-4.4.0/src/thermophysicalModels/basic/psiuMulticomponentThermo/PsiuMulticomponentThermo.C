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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "PsiuMulticomponentThermo.H"
#include "mixtures/inhomogeneousMixture/inhomogeneousThermoMixture.H"
#include "mixtures/homogeneousMixture/homogeneousThermoMixture.H"
#include "mixtures/veryInhomogeneousMixture/veryInhomogeneousThermoMixture.H"
#include "mixtures/egrMixture/egrThermoMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::material
(
    const word& name,
    label mixtureType,
    label patchi
) const
{
    // Internal fields
    if (patchi == -1)
    {
        if (mixtureModel_ == inhomogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template
                lookupObject<inhomogeneousThermoMixture>
                (
                    name + "Model"
                ).inhomogeneousReactantMixtureInternal();
            }
            else
            {
                return matDb().template
                lookupObject<inhomogeneousThermoMixture>
                (
                    name + "Model"
                ).inhomogeneousProductMixtureInternal();
            }
        }
        else if (mixtureModel_ == veryInhomogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template
                lookupObject<veryInhomogeneousThermoMixture>
                (
                    name + "Model"
                ).veryInhomogeneousReactantMixtureInternal();
            }
            else
            {
                return matDb().template
                lookupObject<veryInhomogeneousThermoMixture>
                (
                    name + "Model"
                ).veryInhomogeneousProductMixtureInternal();
            }
        }
        else if (mixtureModel_ == homogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template lookupObject<homogeneousThermoMixture>
                (
                    name + "Model"
                ).homogeneousReactantMixtureInternal();
            }
            else
            {
                return matDb().template lookupObject<homogeneousThermoMixture>
                (
                    name + "Model"
                ).homogeneousProductMixtureInternal();
            }
        }
        else
        {
            if (mixtureType == reactants)
            {
                return matDb().template lookupObject<egrThermoMixture>
                (
                    name + "Model"
                ).egrReactantMixtureInternal();
            }
            else
            {
                return matDb().template lookupObject<egrThermoMixture>
                (
                    name + "Model"
                ).egrProductMixtureInternal();
            }
        }
    }
    else
    {
        if (mixtureModel_ == inhomogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template lookupObject<inhomogeneousThermoMixture>
                (
                    name + "Model"
                ).inhomogeneousReactantMixturePatch(patchi);
            }
            else
            {
                return matDb().template lookupObject<inhomogeneousThermoMixture>
                (
                    name + "Model"
                ).inhomogeneousProductMixturePatch(patchi);
            }
        }
        else if (mixtureModel_ == veryInhomogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template
                lookupObject<veryInhomogeneousThermoMixture>
                (
                    name + "Model"
                ).veryInhomogeneousReactantMixturePatch(patchi);
            }
            else
            {
                return matDb().template
                lookupObject<veryInhomogeneousThermoMixture>
                (
                    name + "Model"
                ).veryInhomogeneousProductMixturePatch(patchi);
            }
        }
        else if (mixtureModel_ == homogeneousThermoMixture::typeName)
        {
            if (mixtureType == reactants)
            {
                return matDb().template lookupObject<homogeneousThermoMixture>
                (
                    name + "Model"
                ).homogeneousReactantMixturePatch(patchi);
            }
            else
            {
                return matDb().template lookupObject<homogeneousThermoMixture>
                (
                    name + "Model"
                ).homogeneousProductMixturePatch(patchi);
            }
        }
        else
        {
            if (mixtureType == reactants)
            {
                return matDb().template lookupObject<egrThermoMixture>
                (
                    name + "Model"
                ).egrReactantMixturePatch(patchi);
            }
            else
            {
                return matDb().template lookupObject<egrThermoMixture>
                (
                    name + "Model"
                ).egrProductMixturePatch(patchi);
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::internalThe
(
    const scalarField& he,
    label mixtureType
) const
{
    volScalarField& heRef = const_cast<volScalarField&>(this->he_);
    const scalarField heBack(this->he_.primitiveField());
    heRef.primitiveFieldRef() = scalarField(he);

    scalarField& Test = this->T_.primitiveFieldRef();

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + this->TRefValue())*1e-4);
    const word CpvName
    (
        (this->he_.name() == "ha" || this->he_.name() == "hs")
      ? CpModel::typeName
      : CvModel::typeName
    );

    int iter = 0;
    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (material(this->he_.name(), mixtureType) - he)
            /material(CpvName, mixtureType);
        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    // Restore temperature in registry
    Test = Told;
    heRef.primitiveFieldRef() = heBack;

    return tTnew;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::boundaryThe
(
    const scalarField& he,
    label mixtureType,
    label patchi
) const
{
    volScalarField& heRef = const_cast<volScalarField&>(this->he_);
    const scalarField heBack(this->he_.boundaryField()[patchi]);
    heRef.boundaryFieldRef()[patchi].forceAssign(scalarField(he));

    scalarField& Test = this->T_.boundaryFieldRef()[patchi];

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + this->TRefValue())*1e-4);
    const word CpvName
    (
        (this->he_.name() == "ha" || this->he_.name() == "hs")
      ? CpModel::typeName
      : CvModel::typeName
    );

    int iter = 0;
    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (material(this->he_.name(), mixtureType, patchi) - he)
            /material(CpvName, mixtureType, patchi);
        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    // Restore temperature in registry
    Test = Told;
    heRef.boundaryFieldRef()[patchi].forceAssign(heBack);

    return tTnew;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const materialTables& mat = this->materials_;
    this->T_.primitiveFieldRef() =
        mat(TModel::typeName, this->phaseName_).primitiveField();

    // Update he and T boundary fields
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he_.boundaryFieldRef();
    forAll(TBf, patchi)
    {
        if (TBf[patchi].fixesValue())
        {
            heBf[patchi].forceAssign
            (
                mat(heModel::typeName, this->phaseName_).boundaryField()[patchi]
            );
        }
        else
        {
            TBf[patchi].forceAssign
            (
                mat(TModel::typeName, this->phaseName_).boundaryField()[patchi]
            );
        }
    }
    this->Cp_.forceAssign(mat(CpModel::typeName, this->phaseName_)());
    this->Cv_.forceAssign(mat(CvModel::typeName, this->phaseName_)());
    this->rho_.forceAssign(mat(rhoModel::typeName, this->phaseName_)());
    this->psi_.forceAssign(mat(psiModel::typeName, this->phaseName_)());
    this->mu_.forceAssign(mat(muModel::typeName, this->phaseName_)());
    this->kappa_.forceAssign(mat(kappaModel::typeName, this->phaseName_)());
    this->Tu_.primitiveFieldRef() =
        internalThe(this->heu_.primitiveField(), reactants);

    forAll(this->T_.boundaryField(), patchi)
    {
        const bool fixesValue =
            this->T_.boundaryFieldRef()[patchi].fixesValue();

        if (!fixesValue)
        {
            this->Tu_.boundaryFieldRef()[patchi].forceAssign
            (
                boundaryThe
                (
                    this->heu_.boundaryField()[patchi],
                    reactants,
                    patchi
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::
PsiuMulticomponentThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo<BasicPsiThermo, MixtureType>(obr, phaseName),
    Tu_
    (
        IOobject
        (
            "Tu",
            obr.time().timeName(),
            obr,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(obr)
    ),
    heu_
    (
        IOobject
        (
            this->he_.name() + 'u',
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimensionSet(0, 2, -2, 0, 0),
        this->heuBoundaryTypes()
    ),
    mixtureModel_(this->properties().template lookup<word>("mixture"))
{
    const scalarField TBack(this->T_.primitiveField());
    this->T_.primitiveFieldRef() = Tu_.primitiveField();
    this->heu_.primitiveFieldRef() = material(this->he_.name(), reactants);
    this->T_.primitiveFieldRef() = TBack;

    volScalarField::Boundary& heuBf = heu_.boundaryFieldRef();
    forAll(heuBf, patchi)
    {
        const scalarField TPatchBack(this->T_.boundaryField()[patchi]);
        this->T_.boundaryFieldRef()[patchi].forceAssign
        (
            this->Tu_.boundaryField()[patchi]
        );
        this->heu_.boundaryFieldRef()[patchi].forceAssign
        (
            material(this->he_.name(), reactants, patchi)
        );
        this->T_.boundaryFieldRef()[patchi].forceAssign(TPatchBack);
    }

    this->heuBoundaryCorrection(this->heu_);

    calculate();
    this->psi_.oldTime();   // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::
~PsiuMulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& p,
    const scalarField& Tu,
    const labelList& cells
) const
{
    tmp<scalarField> theu(new scalarField(cells.size()));
    scalarField& heu = theu.ref();

    const scalarField pOld(this->p_, cells);
    const scalarField TOld(this->T_, cells);
    UIndirectList<scalar>(this->p_, cells) = p;
    UIndirectList<scalar>(this->T_, cells) = Tu;
    // All internal field
    const scalarField internalHe(material(this->he_.name(), reactants));
    forAll(cells, celli)
    {
        heu[celli] = internalHe[cells[celli]];
    }

    // Reset T/p back to the original
    UIndirectList<scalar>(this->T_, cells) = TOld;
    UIndirectList<scalar>(this->p_, cells) = pOld;

    return theu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& Tu,
    const label patchi
) const
{
    const scalarField TBack(this->T_.boundaryField()[patchi]);
    this->T_.boundaryFieldRef()[patchi].forceAssign(Tu);
    tmp<scalarField> theu(material(this->he_.name(), reactants, patchi));
    this->T_.boundaryFieldRef()[patchi].forceAssign(TBack);
    return theu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::Tb() const
{
    tmp<volScalarField> tTb
    (
        volScalarField::New
        (
            "Tb",
            this->T_.db(),
            this->mesh(this->T_.db()),
            dimensionedScalar(dimTemperature, 0)
        )
    );
    volScalarField& tb = tTb.ref();
    tb.primitiveFieldRef() =
        internalThe(this->he_.primitiveField(), products);

    forAll(tb.boundaryField(), patchi)
    {
        tb.boundaryFieldRef()[patchi] =
            boundaryThe(this->he_.boundaryField()[patchi], products, patchi);
    }
    return tTb;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::psiu() const
{
    tmp<volScalarField> tPsiu
    (
        volScalarField::New
        (
            "psiu",
            this->T_.db(),
            this->mesh(this->T_.db()),
            dimensionedScalar(this->psi_.dimensions(), 0)
        )
    );
    volScalarField& psiu = tPsiu.ref();

    const scalarField TOld(this->T_.primitiveField());
    this->T_.primitiveFieldRef() = this->Tu_.primitiveField();
    psiu.primitiveFieldRef() = material(psiModel::typeName, reactants);
    this->T_.primitiveFieldRef() = TOld;

    forAll(psiu.boundaryField(), patchi)
    {
        const scalarField patchTOld(this->T_.boundaryField()[patchi]);
        this->T_.boundaryFieldRef()[patchi].forceAssign
        (
            this->Tu_.boundaryField()[patchi]
        );
        psiu.boundaryFieldRef()[patchi] =
            material(psiModel::typeName, reactants, patchi);
        this->T_.boundaryFieldRef()[patchi].forceAssign(patchTOld);
    }
    return tPsiu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::psib() const
{
    const volScalarField Tb(this->Tb());
    tmp<volScalarField> tPsib
    (
        volScalarField::New
        (
            "psib",
            this->T_.db(),
            this->mesh(this->T_.db()),
            dimensionedScalar(this->psi_.dimensions(), 0)
        )
    );
    volScalarField& psib = tPsib.ref();

    const scalarField TOld(this->T_.primitiveField());
    this->T_.primitiveFieldRef() = Tb.primitiveField();
    psib.primitiveFieldRef() = material(psiModel::typeName, products);
    this->T_.primitiveFieldRef() = TOld;

    forAll(psib.boundaryField(), patchi)
    {
        const scalarField patchTOld(this->T_.boundaryField()[patchi]);
        this->T_.boundaryFieldRef()[patchi].forceAssign
        (
            Tb.boundaryField()[patchi]
        );
        psib.boundaryFieldRef()[patchi] =
            material(psiModel::typeName, products, patchi);
        this->T_.boundaryFieldRef()[patchi].forceAssign(patchTOld);
    }
    return tPsib;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::muu() const
{
    tmp<volScalarField> tmmu
    (
        volScalarField::New
        (
            "muu",
            this->T_.db(),
            this->mesh(this->T_.db()),
            dimensionedScalar(dimDynamicViscosity, 0)
        )
    );
    volScalarField& mmu = tmmu.ref();

    const scalarField TOld(this->T_.primitiveField());
    this->T_.primitiveFieldRef() = Tu_.primitiveField();
    mmu.primitiveFieldRef() = material(muModel::typeName, products);
    this->T_.primitiveFieldRef() = TOld;

    forAll(mmu.boundaryField(), patchi)
    {
        const scalarField patchTOld(this->T_.boundaryField()[patchi]);
        this->T_.boundaryFieldRef()[patchi].forceAssign
        (
            Tu_.boundaryField()[patchi]
        );
        mmu.boundaryFieldRef()[patchi] =
            material(muModel::typeName, reactants, patchi);
        this->T_.boundaryFieldRef()[patchi].forceAssign(patchTOld);
    }
    return tmmu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>::mub() const
{
    const volScalarField Tb(this->Tb());
    tmp<volScalarField> tmub
    (
        volScalarField::New
        (
            "mub",
            this->T_.db(),
            this->mesh(this->T_.db()),
            dimensionedScalar(dimDynamicViscosity, 0)
        )
    );
    volScalarField& mmb = tmub.ref();

    const scalarField TOld(this->T_.primitiveField());
    this->T_.primitiveFieldRef() = Tb.primitiveField();
    mmb.primitiveFieldRef() = material(muModel::typeName, products);
    this->T_.primitiveFieldRef() = TOld;

    forAll(mmb.boundaryField(), patchi)
    {
        const scalarField patchTOld(this->T_.boundaryField()[patchi]);
        this->T_.boundaryFieldRef()[patchi].forceAssign
        (
            Tb.boundaryField()[patchi]
        );
        mmb.boundaryFieldRef()[patchi] =
            material(muModel::typeName, products, patchi);
        this->T_.boundaryFieldRef()[patchi].forceAssign(patchTOld);
    }
    return tmub;
}


// ************************************************************************* //

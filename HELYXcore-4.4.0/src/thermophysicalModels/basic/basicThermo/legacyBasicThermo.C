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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "basicThermo/legacyBasicThermo.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
template
<
    class CellMixture,
    class PatchFaceMixture,
    class Method,
    class ... Args
>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    CellMixture cellMixture,
    PatchFaceMixture patchFaceMixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->phaseName_),
            this->T_.db(),
            this->T_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->T_, celli)
    {
        psi[celli] = ((this->*cellMixture)(celli).*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->T_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                ((this->*patchFaceMixture)(patchi, facei).*psiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                );
        }
    }

    return tPsi;
}


template<class BasicThermo, class MixtureType>
template<class PatchFaceMixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::patchFieldProperty
(
    PatchFaceMixture patchFaceMixture,
    Method psiMethod,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->T_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        psi[facei] =
            ((this->*patchFaceMixture)(patchi, facei).*psiMethod)
            (
                args[facei] ...
            );
    }

    return tPsi;
}

template<class BasicThermo, class MixtureType>
void Foam::legacyBasicThermo<BasicThermo, MixtureType>::init
(
    const volScalarField& p,
    const volScalarField& T,
    volScalarField& he
)
{
    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p.primitiveField();
    const scalarField& TCells = T.primitiveField();

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellThermoMixture(celli).he
            (
                pCells[celli] + this->pRefValue(),
                TCells[celli]
            );
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        heBf[patchi].forceAssign(this->he(T.boundaryField()[patchi], patchi));
    }

    this->heBoundaryCorrection(he);

    // Note: T does not have oldTime
    if (p.nOldTimes() > 0)
    {
        init(p.oldTime(), T.oldTime(), he.oldTime());
    }
}


template<class BasicThermo, class MixtureType>
void Foam::legacyBasicThermo<BasicThermo, MixtureType>::heBoundaryCorrection
(
    fvPatchScalarField& pf
)
{
    if (isA<gradientEnergyFvPatchScalarField>(pf))
    {
        refCast<gradientEnergyFvPatchScalarField>(pf).gradient() =
            pf.fvPatchField::snGrad();
    }
    else if (isA<mixedEnergyFvPatchScalarField>(pf))
    {
        refCast<mixedEnergyFvPatchScalarField>(pf).refGrad() =
            pf.fvPatchField::snGrad();
    }
    else if (isA<blendedEnergyFvPatchScalarField>(pf))
    {
        blendedEnergyFvPatchScalarField& bepf =
            refCast<blendedEnergyFvPatchScalarField>(pf);
        heBoundaryCorrection(bepf.boundaryOne());
        heBoundaryCorrection(bepf.boundaryTwo());
    }
}


template<class BasicThermo, class MixtureType>
void Foam::legacyBasicThermo<BasicThermo, MixtureType>::heBoundaryCorrection
(
    volScalarField& h
)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        heBoundaryCorrection(hBf[patchi]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::legacyBasicThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo(obr, phaseName),
    MixtureType(*this, obr, phaseName),
    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    Cp_
    (
        IOobject
        (
            BasicThermo::phasePropertyName("thermo-Cp", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            BasicThermo::phasePropertyName("thermo-Cv", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    init(this->p_, this->T_, he_);
}


template<class BasicThermo, class MixtureType>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::legacyBasicThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(obr, dict, phaseName),
    MixtureType(*this, obr, phaseName),
    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    Cp_
    (
        IOobject
        (
            BasicThermo::phasePropertyName("thermo-Cp", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            BasicThermo::phasePropertyName("thermo-Cv", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    init(this->p_, this->T_, he_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::~legacyBasicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::legacyBasicThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const tmp<volScalarField> pAbsolute
    (
        volScalarField::New("pAbs", p + this->pRef())
    );
    return volScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::he,
        pAbsolute(),
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::legacyBasicThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, celli)
    {
        he[celli] =
            this->cellThermoMixture(cells[celli]).he
            (
                p[celli] + this->pRefValue(),
                T[celli]
            );
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::legacyBasicThermo<BasicThermo, MixtureType>::he
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::he,
        patchi,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::hs() const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::hs,
        this->pAbs()(),
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::hs,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::hs
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> ths(new scalarField(T.size()));
    forAll(T, celli)
    {
        ths.ref()[celli] =
            this->cellThermoMixture(cells[celli]).hs
            (
                p[celli] + this->pRefValue(),
                T[celli]
            );
    }
    return ths;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::hs,
        patchi,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::ha() const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::ha,
        this->pAbs()(),
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const tmp<volScalarField> pAbsolute
    (
        volScalarField::New("pAbs", p + this->pRef())
    );
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::ha,
        pAbsolute(),
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::ha,
        patchi,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::hf() const
{
    return volScalarFieldProperty
    (
        "hf",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::hf
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cp,
        patchi,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cv,
        patchi,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return Cp(T, patchi);
    }
    else
    {
        return Cv(T, patchi);
    }
}


template<class BasicThermo, class MixtureType>
const Foam::volScalarField&
Foam::legacyBasicThermo<BasicThermo, MixtureType>::Cpv() const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return Cp_;
    }
    else
    {
        return Cv_;
    }
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::The
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    const scalarField p(this->p().boundaryField()[patchi] + this->pRefValue());
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::The,
        patchi,
        h,
        p,
        T0
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::The
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return volScalarFieldProperty
    (
        "T",
        dimTemperature,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::The,
        h,
        p,
        T0
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New("kappaEff", this->kappa_ + Cp_*alphat);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->kappa_.boundaryField()[patchi]
      + Cp(this->T_.boundaryField()[patchi], patchi)*alphat;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return volScalarField::New("alphaEff", this->kappa_/Cp_ + alphat);
    }
    else
    {
        return volScalarField::New
        (
            "alphaEff",
            (this->kappa_ + Cp_*alphat)/Cv_
        );
    }
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return
            this->kappa_.boundaryField()[patchi]/Cp_.boundaryField()[patchi]
          + alphat;
    }
    else
    {
        return
            (
                this->kappa_.boundaryField()[patchi]
              + Cp_.boundaryField()[patchi]*alphat
            )/Cv_.boundaryField()[patchi];
    }
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::legacyBasicThermo<BasicThermo, MixtureType>::W() const
{
    return volScalarFieldProperty
    (
        "W",
        dimMass/dimMoles,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::W
    );
}


template<class BasicThermo, class MixtureType>
bool Foam::legacyBasicThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(this->mesh(BasicThermo::db()), *this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

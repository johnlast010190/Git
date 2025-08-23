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
    (c) 2013-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "twoPhaseMixtureThermo.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::twoPhaseMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rhoThermo::composite(U.mesh(), word::null),
    twoPhaseMixture(U.mesh(), *this),
    interfaceProperties(alpha1(), U, *this),
    thermo1_(nullptr),
    thermo2_(nullptr),
    Cp_
    (
        IOobject
        (
            "thermo-Cp",
            U.db().time().timeName(),
            U.db()
        ),
        this->mesh(U.db()),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            "thermo-Cv",
            U.db().time().timeName(),
            U.db()
        ),
        this->mesh(U.db()),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    {
        volScalarField T1(IOobject::groupName("T", phase1Name()), T_);
        T1.write();
    }

    {
        volScalarField T2(IOobject::groupName("T", phase2Name()), T_);
        T2.write();
    }

    thermo1_ = rhoThermo::New(U.mesh(), phase1Name());
    thermo2_ = rhoThermo::New(U.mesh(), phase2Name());

    thermo1_->validate(phase1Name(), "e");
    thermo2_->validate(phase2Name(), "e");

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::~twoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseMixtureThermo::correctThermo()
{
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}


void Foam::twoPhaseMixtureThermo::correct()
{
    rho_ = alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho();
    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    kappa_ = alpha1()*thermo1_->kappa() + alpha2()*thermo2_->kappa();
    Cp_ = alpha1()*thermo1_->Cp() + alpha2()*thermo2_->Cp();
    Cv_ = alpha1()*thermo1_->Cv() + alpha2()*thermo2_->Cv();

    interfaceProperties::correct();
}


Foam::word Foam::twoPhaseMixtureThermo::thermoName() const
{
    return thermo1_->thermoName() + ',' + thermo2_->thermoName();
}


bool Foam::twoPhaseMixtureThermo::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}


bool Foam::twoPhaseMixtureThermo::isochoric() const
{
    return thermo1_->isochoric() && thermo2_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->he(p, T) + alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->he(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->he(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hs() const
{
    return alpha1()*thermo1_->hs() + alpha2()*thermo2_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->hs(p, T) + alpha2()*thermo2_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::hs
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->hs(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->hs(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->hs(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->hs(T, patchi);
}



Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::ha() const
{
    return alpha1()*thermo1_->ha() + alpha2()*thermo2_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->ha(p, T) + alpha2()*thermo2_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->ha(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hf() const
{
    return alpha1()*thermo1_->hf() + alpha2()*thermo2_->hf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::The
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::The
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
}


const Foam::volScalarField& Foam::twoPhaseMixtureThermo::Cp() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cp(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cp(T, patchi);
}


const Foam::volScalarField& Foam::twoPhaseMixtureThermo::Cv() const
{
    return Cv_;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cv(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cv(T, patchi);
}



const Foam::volScalarField& Foam::twoPhaseMixtureThermo::Cpv() const
{
    if (thermo1_->he().member() == "h" || thermo1_->he().member() == "ha")
    {
        return Cp_;
    }
    else
    {
        return Cv_;
    }
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cpv(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::nu() const
{
    return mu()/(alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho());
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::nu
(
    const label patchi
) const
{
    return
        mu().boundaryField()[patchi]
       /(
            alpha1().boundaryField()[patchi]*thermo1_->rho(patchi)
          + alpha2().boundaryField()[patchi]*thermo2_->rho(patchi)
        );
}


const Foam::volScalarField& Foam::twoPhaseMixtureThermo::kappa() const
{
    return kappa_;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return kappa_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->kappaEff(alphat)
      + alpha2()*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->alphaEff(alphat)
      + alpha2()*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi)
    ;
}


bool Foam::twoPhaseMixtureThermo::read()
{
    if (rhoThermo::composite::read())
    {
        return interfaceProperties::read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

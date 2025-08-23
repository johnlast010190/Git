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
    (c) 2015-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "AnisothermalPhaseModel.H"
#include "../../eulerianPhaseSystem/eulerianPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::AnisothermalPhaseModel
(
    const eulerianPhaseSystem& fluid,
    const word& phaseName,
    const label index,
    const dictionary& coeffs
)
:
    BasePhaseModel(fluid, phaseName, index, coeffs),
    K_
    (
        IOobject
        (
            IOobject::groupName("K", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(sqr(dimVelocity), 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::~AnisothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();
    K_ = 0.5*magSqr(this->U());
}


template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    // Main thermo is already corrected in eulerianPhaseSystem
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->properties().lookupOrDefault
        (
            "pressureWorkAlphaLimit",
            0.0
        );

    if (pressureWorkAlphaLimit > 0)
    {
        return
        (
            max(alpha - pressureWorkAlphaLimit, scalar(0))
           /max(alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)
        )*pressureWork;
    }
    else
    {
        return pressureWork;
    }
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;
    const volVectorField& U = this->U();
    tmp<surfaceScalarField> talphaPhiv(this->alphaPhiv());
    const surfaceScalarField& alphaPhiv = talphaPhiv();
    const surfaceScalarField& alphaPhi = this->alphaPhi();

    const tmp<volScalarField> tcontErr(this->continuityError());
    const volScalarField& contErr(tcontErr());

    const volScalarField alphaEff(this->turbulence().alphaEff());

    volScalarField& he = this->thermo_->he();

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, this->rho(), he)
      + fvm::div(alphaPhi, he)
      - fvm::Sp(contErr, he)

      + fvc::ddt(alpha, this->rho(), K_) + fvc::div(alphaPhi, K_)
      - contErr*K_

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(alphaEff),
            he
        )
     ==
        this->Qdot()
    );

    // Add the appropriate pressure-work term
    if (he.name() == this->thermo_->phasePropertyName("e"))
    {
        tEEqn.ref() += filterPressureWork
        (
            fvc::div(fvc::absolute(alphaPhiv, alpha, U), this->thermo().p())
          + this->thermo().p()*fvc::ddt(alpha)
        );
    }
    else if (this->thermo_->dpdt())
    {
        tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
    }

    return tEEqn;
}


template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::compressible() const
{
    return !this->thermo().incompressible();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::AnisothermalPhaseModel<BasePhaseModel>::K() const
{
    return K_;
}


// ************************************************************************* //

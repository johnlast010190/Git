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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "MovingPhaseModel.H"
#include "../../eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/slip/slipFvPatchFields.H"
#include "fields/fvPatchFields/derived/partialSlip/partialSlipFvPatchFields.H"

#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phiv(const volVectorField& U) const
{
    word phiName(IOobject::groupName("phiv", this->name()));

    IOobject phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U),
                phiTypes
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::MovingPhaseModel
(
    const eulerianPhaseSystem& fluid,
    const word& phaseName,
    const label index,
    const dictionary& coeffs
)
:
    BasePhaseModel(fluid, phaseName, index, coeffs),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    phiv_(phiv(U_)),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    DUDt_
    (
        IOobject
        (
            IOobject::groupName("DUDt", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedVector(dimAcceleration, Zero)
    ),
    divU_(nullptr),
    turbulence_
    (
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            this->thermo().rho(),
            U_,
            alphaPhi_,
            phiv_,
            *this
        )
    ),
    continuityError_
    (
        IOobject
        (
            IOobject::groupName("continuityError", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    )
{
    alphaPhi_.setOriented();

    phiv_.writeOpt() = IOobject::AUTO_WRITE;
    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();

    volScalarField& rho = this->thermo().rho();

    continuityError_ =
        fvc::ddt(*this, rho) + fvc::div(alphaPhi_)
      - (this->fluid().fvOptions()(*this, rho)&rho);
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();

    DUDt_ = fvc::ddt(U_) + fvc::div(phiv_, U_) - fvc::div(phiv_)*U_;
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctTurbulence()
{
    BasePhaseModel::correctTurbulence();

    turbulence_->correct();
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctEnergyTransport()
{
    BasePhaseModel::correctEnergyTransport();
    turbulence_->correctEnergyTransport();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UEqn()
{
    return
    (
        fvm::ddt(*this, this->thermo().rho(), U_)
      + fvm::div(alphaPhi_, U_)
      - fvm::Sp(continuityError_, U_)
      + this->fluid().fvOptions().MRFDDt(this->volFrac()*this->thermo().rho(), U_)
      + turbulence_->divDevRhoReff(U_)
    );
}


template<class BasePhaseModel>
const Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::DbyA() const
{
    if (DbyA_.valid())
    {
        return DbyA_;
    }
    else
    {
        return surfaceScalarField::null();
    }
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::DbyA
(
    const tmp<surfaceScalarField>& DbyA
)
{
    DbyA_ = DbyA;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::U() const
{
    return U_;
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::MovingPhaseModel<BasePhaseModel>::U()
{
    return U_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::DUDt() const
{
    return DUDt_;
}


template<class BasePhaseModel>
const Foam::tmp<Foam::volScalarField>&
Foam::MovingPhaseModel<BasePhaseModel>::divU() const
{
    return divU_;
}


template<class BasePhaseModel>
void
Foam::MovingPhaseModel<BasePhaseModel>::divU
(
    const tmp<volScalarField>& divU
)
{
    divU_ = divU;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityError() const
{
    return continuityError_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phiv() const
{
    return phiv_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::phiv()
{
    return phiv_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhiv() const
{
    return
        phiv_.db().template lookupObjectRef<surfaceScalarField>
        (
            IOobject::groupName("alphaPhiv", this->name())
        );
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi()
{
    return alphaPhi_;
}


template<class BasePhaseModel>
const Foam::phaseCompressibleTurbulenceModel&
Foam::MovingPhaseModel<BasePhaseModel>::turbulence() const
{
    return turbulence_;
}


// ************************************************************************* //

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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2023 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "eulerianPhaseModel.H"
#include "../../eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "../../eulerianDiameterModels/eulerianDiameterModel/eulerianDiameterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianPhaseModel, 0);
    defineRunTimeSelectionTable(eulerianPhaseModel, eulerianPhaseSystem);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianPhaseModel::eulerianPhaseModel
(
    const eulerianPhaseSystem& fluid,
    const word& phaseName,
    const label index,
    const dictionary& coeffs
)
:
    fluid_(fluid),
    name_(phaseName),
    index_(index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffs.lookup("residualAlpha")
    ),
    alphaMax_(coeffs.lookupOrDefault("alphaMax", 1.0))
{
    diameterModel_ = eulerianDiameterModel::New(coeffs, *this);
}


Foam::autoPtr<Foam::eulerianPhaseModel> Foam::eulerianPhaseModel::clone() const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianPhaseModel::~eulerianPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::eulerianPhaseModel>
Foam::eulerianPhaseModel::iNew::operator()(Istream& is) const
{
    indexCounter_++;
    word phaseName(is);
    return autoPtr<eulerianPhaseModel>
    (
        eulerianPhaseModel::New
        (
            fluid_,
            phaseName,
            indexCounter_,
            fluid_.thermo().thermos(phaseName).phaseDict()
        )
    );
}


Foam::eulerianPhaseModel::operator const volScalarField&() const
{
    return fluid_.alphas()[index_];
}


Foam::eulerianPhaseModel::operator volScalarField&()
{
    return const_cast<eulerianPhaseSystem&>(fluid_).alphas()[index_];
}


const Foam::word& Foam::eulerianPhaseModel::name() const
{
    return name_;
}


const Foam::word& Foam::eulerianPhaseModel::keyword() const
{
    return name_;
}


Foam::label Foam::eulerianPhaseModel::index() const
{
    return index_;
}


const Foam::eulerianPhaseSystem& Foam::eulerianPhaseModel::fluid() const
{
    return fluid_;
}


const Foam::dimensionedScalar& Foam::eulerianPhaseModel::residualAlpha() const
{
    return residualAlpha_;
}


Foam::scalar Foam::eulerianPhaseModel::alphaMax() const
{
    return alphaMax_;
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhaseModel::d() const
{
    return diameterModel_().d();
}


void Foam::eulerianPhaseModel::correct()
{
    diameterModel_->correct();
}


void Foam::eulerianPhaseModel::correctKinematics()
{}


void Foam::eulerianPhaseModel::correctThermo()
{}


void Foam::eulerianPhaseModel::correctTurbulence()
{}


void Foam::eulerianPhaseModel::correctEnergyTransport()
{}


bool Foam::eulerianPhaseModel::read()
{
    return
        diameterModel_->read
        (
            fluid_.thermo().thermos(name_).phaseDict().subDict
            (
                this->typeName + "Coeffs"
            )
        );

}


bool Foam::eulerianPhaseModel::compressible() const
{
    return false;
}


void Foam::eulerianPhaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhivBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = volFrac().boundaryField();

    tmp<surfaceScalarField> tphiv(this->phiv());
    const surfaceScalarField& phiv = tphiv();
    const surfaceScalarField::Boundary& phivBf = phiv.boundaryField();

    forAll(alphaPhivBf, patchi)
    {
        fvsPatchScalarField& alphaPhivp = alphaPhivBf[patchi];

        if (!alphaPhivp.coupled())
        {
            alphaPhivp = phivBf[patchi]*alphaBf[patchi];
        }
    }
}


const Foam::tmp<Foam::volScalarField>& Foam::eulerianPhaseModel::divU() const
{
    NotImplemented;
}


void Foam::eulerianPhaseModel::divU(const tmp<volScalarField>& divU)
{
    WarningInFunction
        << "Attempt to set the dilatation rate of an incompressible phase"
        << endl;
}


const Foam::volScalarField& Foam::eulerianPhaseModel::K() const
{
    NotImplemented;
}


const Foam::surfaceScalarField& Foam::eulerianPhaseModel::DbyA() const
{
    return surfaceScalarField::null();
}


void Foam::eulerianPhaseModel::DbyA(const tmp<surfaceScalarField>& DbyA)
{
    WarningInFunction
        << "Attempt to set the dilatation rate of an incompressible phase"
        << endl;
}


// ************************************************************************* //

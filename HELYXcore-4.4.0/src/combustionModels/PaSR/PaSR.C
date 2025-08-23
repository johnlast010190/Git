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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "PaSR/PaSR.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(PaSR, 0);
    addToRunTimeSelectionTable(combustionModel, PaSR, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::PaSR
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, turb, combustionProperties),
    Cmix_(this->coeffs().template lookup<scalar>("Cmix")),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":kappa"),
            this->thermo().T().db().time().timeName(),
            this->thermo().T().db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::PaSR::correct()
{
    laminar::correct();

    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    const scalarField& epsilon = tepsilon();

    tmp<volScalarField> tmuEff(this->turbulence().muEff());
    const scalarField& muEff = tmuEff();

    tmp<volScalarField> ttc(this->chemistryPtr_->tc());
    const scalarField& tc = ttc();

    tmp<volScalarField> trho(this->rho());
    const scalarField& rho = trho();

    forAll(epsilon, i)
    {
        const scalar tk =
            Cmix_*sqrt(max(muEff[i]/rho[i]/(epsilon[i] + SMALL), 0));
        if (tk > SMALL)
        {
            kappa_[i] = tc[i]/(tc[i] + tk);
        }
        else
        {
            kappa_[i] = 1.0;
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR::R(volScalarField& Y) const
{
    return kappa_*laminar::R(Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + "-Qdot"),
        kappa_*laminar::Qdot()
    );
}


bool Foam::combustionModels::PaSR::read()
{
    if (laminar::read())
    {
        Cmix_ = this->coeffs().template lookup<scalar>("Cmix");
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

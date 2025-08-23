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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2025 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "phaseChangeModel.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, multiphase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const dictionary& dict,
    const phasePairKey& pair,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    pair_(pair),
    thermo_(refCast<multiphaseThermo>(multiphaseThermo::lookupOrCreate(mesh_))),
    phaseIdx1_(thermo_.phases().find(pair.first())),
    phaseIdx2_(thermo_.phases().find(pair.second())),
    pSat_("pSat", dimPressure, 1.0),
    pSatTable_(Function1<scalar>::New("pSat", dict))
{
    pSat_.value() = pSatTable_->value(mesh_.time().value());
    if (phaseIdx1_ < 0)
    {
        FatalErrorInFunction
            << "Phase " << pair.first() << " was not found."
            << nl << exit(FatalError);
    }
    if (phaseIdx2_ < 0)
    {
        FatalErrorInFunction
            << "Phase " << pair.second() << " was not found."
            << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::phaseChangeModel::alpha1() const
{
    return thermo_.alphas()[phaseIdx1_];
}

const Foam::volScalarField& Foam::phaseChangeModel::alpha2() const
{
    return thermo_.alphas()[phaseIdx2_];
}

Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::rho1() const
{
    return thermo_.thermos()[phaseIdx1_].rho();
}

Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::rho2() const
{
    return thermo_.thermos()[phaseIdx2_].rho();
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModel::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1()*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModel::vDotP() const
{
    volScalarField pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


void Foam::phaseChangeModel::correct()
{
    pSat_.value() = pSatTable_->value(mesh_.time().value());
}

// ************************************************************************* //

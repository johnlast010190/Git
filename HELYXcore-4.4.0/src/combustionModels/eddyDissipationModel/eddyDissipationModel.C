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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "eddyDissipationModel/eddyDissipationModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(eddyDissipationModel, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        eddyDissipationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationModel::eddyDissipationModel
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion(modelType, thermo, turb, combustionProperties),
    C_(this->coeffs().template lookup<scalar>("CEDC")),
    Cd_(this->coeffs().template lookup<scalar>("CDiff")),
    Cstiff_(this->coeffs().template lookup<scalar>("CStiff")),
    PV_
    (
        IOobject
        (
            "PV",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationModel::~eddyDissipationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationModel::rtTurb() const
{
    return
        C_*this->turbulence().epsilon()
       /max
        (
            this->turbulence().k(),
            dimensionedScalar(dimVelocity*dimVelocity, SMALL)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationModel::rtDiff() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );

    return
        Cd_*(this->thermo().kappa()/this->thermo().Cp())
       /this->rho()/sqr(lesModel.delta());
}


void Foam::combustionModels::eddyDissipationModel::correct()
{
    // Set the product volume field, needed by alphat BC
    calcPV();

    this->wFuel_.forceAssign
    (
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0)
    );

    this->fresCorrect();
    const label fuelI = this->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    if (this->thermo().composition().contains("O2"))
    {
        const volScalarField& YO2 = this->thermo().composition().Y("O2");
        volScalarField rt(max(rtTurb(), rtDiff()));
        this->wFuel_.forceAssign
        (
            this->rho()
           *min(YFuel, YO2/this->s().value())
           /this->mesh_.time().deltaT()/Cstiff_
           *(1 - exp(-Cstiff_*this->mesh_.time().deltaT()*rt))
        );
    }
}


bool Foam::combustionModels::eddyDissipationModel::read()
{
    if (singleStepCombustion::read())
    {
        C_ = this->coeffs().template lookup<scalar>("C");
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::combustionModels::eddyDissipationModel::calcPV()
{
    // Get species mass fraction
    const label fuelI = this->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const volScalarField& YCO2 = this->thermo().composition().Y("CO2");
    const scalar s = this->s().value();

    // Get Mspecies/Mfuel from reaction equation
    const label CO2i = this->thermo().composition().species()["CO2"];
    const scalar rCO2(this->specieStoichCoeffs()[CO2i]);
    const label H2Oi = this->thermo().composition().species()["CO2"];
    const scalar rH2O(this->specieStoichCoeffs()[H2Oi]);

    PV_ =
        (YCO2*(1.0 + rH2O/rCO2) + SMALL)
       /(YCO2*(1.0 + rH2O/rCO2) + SMALL + min(YFuel, YO2/s)*(1.0 + s));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

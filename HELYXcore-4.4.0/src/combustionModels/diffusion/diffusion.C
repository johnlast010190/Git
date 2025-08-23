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
    (c) 2012-2022 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "diffusion/diffusion.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(diffusion, 0);
    addToRunTimeSelectionTable(combustionModel, diffusion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::diffusion::diffusion
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion(modelType, thermo, turb, combustionProperties),
    C_(this->coeffs().template lookup<scalar>("C")),
    oxidantName_(this->coeffs().template lookupOrDefault<word>("oxidant", "O2"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::diffusion::~diffusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::diffusion::correct()
{
    this->wFuel_.forceAssign
    (
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0)
    );

    this->fresCorrect();

    const label fuelI = this->fuelIndex();

    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];

    if (this->thermo().composition().contains(oxidantName_))
    {
        const volScalarField& YO2 =
            this->thermo().composition().Y(oxidantName_);
        this->wFuel_.forceAssign
        (
            C_*this->turbulence().muEff()
           *mag(fvc::grad(YFuel) & fvc::grad(YO2))
           *pos0(YFuel)*pos0(YO2)
        );
    }
}


bool Foam::combustionModels::diffusion::read()
{
    if (singleStepCombustion::read())
    {
        C_ = this->coeffs().template lookup<scalar>("C");
        this->coeffs().readIfPresent("oxidant", oxidantName_);
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

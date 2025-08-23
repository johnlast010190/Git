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
    (c) 2016 OpenCFD Ltd
    (c) 2024-2025 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "eddyDissipationModelBase/eddyDissipationModelBase.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationModelBase::eddyDissipationModelBase
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion(modelType, thermo, turb, combustionProperties),
    CEDC_(this->coeffs().template lookup<scalar>("CEDC"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationModelBase::~eddyDissipationModelBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationModelBase::rtTurb() const
{
    return
        CEDC_*this->turbulence().epsilon()
       /max
        (
            this->turbulence().k(),
            dimensionedScalar(sqr(dimVelocity), SMALL)
        );
}


void Foam::combustionModels::eddyDissipationModelBase::correct()
{
    this->wFuel_.forceAssign
    (
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0)
    );

    this->fresCorrect();

    const label fuelI = this->fuelIndex();

    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];

    const dimensionedScalar s = this->s();

    if (this->thermo().composition().contains("O2"))
    {
        const volScalarField& YO2 =
            this->thermo().composition().Y("O2");
        this->wFuel_.forceAssign(
              this->rho()
            * min(YFuel, YO2/s.value())
            * timeScale()
        );
    }
    else
    {
        FatalErrorInFunction
            << "You selected a combustion model which requieres O2 mass"
            << " to be present in the mixture"
            << exit(FatalError);
    }
}


bool Foam::combustionModels::eddyDissipationModelBase::read()
{
    if (singleStepCombustion::read())
    {
        CEDC_ = this->coeffs().template lookup<scalar>("CEDC");
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

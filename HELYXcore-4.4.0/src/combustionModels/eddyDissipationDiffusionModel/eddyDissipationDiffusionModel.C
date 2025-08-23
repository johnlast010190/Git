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

#include "eddyDissipationDiffusionModel/eddyDissipationDiffusionModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(eddyDissipationDiffusionModel, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        eddyDissipationDiffusionModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationDiffusionModel::
eddyDissipationDiffusionModel
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    eddyDissipationModelBase(modelType, thermo, turb, combustionProperties),
    Cd_(this->coeffs().template lookup<scalar>("Cd"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationDiffusionModel::
~eddyDissipationDiffusionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationDiffusionModel::timeScale()
{
    return (max(this->rtTurb(), this->rtDiff()));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationDiffusionModel::rtDiff() const
{
    tmp<volScalarField> tdelta
    (
        volScalarField::New
        (
            "tdelta",
            this->mesh(),
            dimensionedScalar(dimLength, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& delta = tdelta.ref();
    delta.ref() = pow(this->mesh().V(), 1.0/3.0);
    delta.correctBoundaryConditions();

    // NOTE: Assume Prt = 1
    return Cd_*this->turbulence().nuEff()/sqr(delta);
}


bool Foam::combustionModels::eddyDissipationDiffusionModel::read()
{
    if (eddyDissipationModelBase::read())
    {
        Cd_ = this->coeffs().template lookup<scalar>("Cd");
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

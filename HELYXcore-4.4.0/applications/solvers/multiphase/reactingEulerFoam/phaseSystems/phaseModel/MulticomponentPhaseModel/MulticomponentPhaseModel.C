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

#include "MulticomponentPhaseModel.H"

#include "phaseSystem/phaseSystem.H"

#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MulticomponentPhaseModel<BasePhaseModel>::MulticomponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    Sc_("Sc", dimless, fluid.subDict(phaseName)),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.mesh().solution().solverDict("Yi")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MulticomponentPhaseModel<BasePhaseModel>::~MulticomponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MulticomponentPhaseModel<BasePhaseModel>::correctThermo()
{
    this->thermo_->composition().normalise();
    BasePhaseModel::correctThermo();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MulticomponentPhaseModel<BasePhaseModel>::YiEqn
(
    volScalarField& Yi
)
{
    const label defalutSpecie = this->thermo_->composition().defaultSpecie();
    if
    (
        (defalutSpecie != -1)
     && (
            (
                Yi.name()
             == IOobject::groupName
                (
                    this->thermo_->composition().species()[defalutSpecie],
                    this->name()
                )
            )
         || (
               !this->thermo_->composition().active
                (
                    this->thermo_->composition().species()[Yi.member()]
                )
            )
        )
    )
    {
        return tmp<fvScalarMatrix>();
    }

    const volScalarField& alpha = *this;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi();
    const tmp<volScalarField> trho(this->rho());
    const volScalarField& rho(trho());

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")
      - fvm::Sp(this->continuityError(), Yi)

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(this->turbulence().nut()*rho/Sc_),
            Yi
        )
     ==
        this->R(Yi)

      + fvc::ddt(residualAlpha_*rho, Yi)
      - fvm::ddt(residualAlpha_*rho, Yi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::Y()
{
    return this->thermo_->composition().Y();
}


// ************************************************************************* //

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

#include "ReactingPhaseModel.H"
#include "../../eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ReactingPhaseModel<BasePhaseModel>::ReactingPhaseModel
(
    const eulerianPhaseSystem& fluid,
    const word& phaseName,
    const label index,
    const dictionary& coeffs
)
:
    BasePhaseModel(fluid, phaseName, index, coeffs),
    reaction_(combustionModel::New(*this->thermo_, this->turbulence_(), phaseName))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ReactingPhaseModel<BasePhaseModel>::~ReactingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::ReactingPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();
    reaction_->correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix> Foam::ReactingPhaseModel<BasePhaseModel>::R
(
    volScalarField& Yi
) const
{
    return reaction_->R(Yi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::ReactingPhaseModel<BasePhaseModel>::Qdot() const
{
    return reaction_->Qdot();
}


// ************************************************************************* //

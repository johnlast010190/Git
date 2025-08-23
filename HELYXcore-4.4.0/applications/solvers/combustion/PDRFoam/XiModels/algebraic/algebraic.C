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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "algebraic.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(algebraic, 0);
    addToRunTimeSelectionTable(XiModel, algebraic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::algebraic
(
    const dictionary& XiProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su,
    const volScalarField& rho,
    const volScalarField& b,
    const surfaceScalarField& phi
)
:
    XiModel(XiProperties, thermo, turbulence, Su, rho, b, phi),
    XiShapeCoeff_(XiModelCoeffs_.lookup<scalar>("XiShapeCoef")),
    XiEqModel_(XiEqModel::New(XiProperties, thermo, turbulence, Su)),
    XiGModel_(XiGModel::New(XiProperties, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::~algebraic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::algebraic::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::algebraic::correct()
{
    volScalarField XiEqEta(XiEqModel_->XiEq());
    volScalarField GEta(XiGModel_->G());

    volScalarField R(GEta*XiEqEta/(XiEqEta - 0.999));

    volScalarField XiEqStar(R/(R - GEta));

    Xi_.forceAssign(1.0 + (1.0 + (2*XiShapeCoeff_)*(0.5 - b_))*(XiEqStar - 1.0));
}


bool Foam::XiModels::algebraic::read(const dictionary& XiProperties)
{
    XiModel::read(XiProperties);
    XiShapeCoeff_ = XiModelCoeffs_.lookup<scalar>("XiShapeCoef");

    return true;
}


// ************************************************************************* //

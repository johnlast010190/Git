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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "basicXiSubG.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(basicSubGrid, 0);
    addToRunTimeSelectionTable(XiGModel, basicSubGrid, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::basicSubGrid::basicSubGrid
(
    const dictionary& XiGProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModel(XiGProperties, thermo, turbulence, Su),
    k1_(XiGModelCoeffs_.lookup<scalar>("k1")),
    XiGModel_(XiGModel::New(XiGModelCoeffs_, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiGModels::basicSubGrid::~basicSubGrid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::basicSubGrid::G() const
{
    const objectRegistry& db = Su_.db();
    const volVectorField& U = db.lookupObject<volVectorField>("U");
    const volScalarField& Nv = db.lookupObject<volScalarField>("Nv");
    const volScalarField& Lobs = db.lookupObject<volScalarField>("Lobs");

    tmp<volScalarField> tGtot = XiGModel_->G();
    volScalarField& Gtot = tGtot.ref();

    const scalarField Cw = pow(Su_.mesh().V(), 2.0/3.0);
    scalarField N(Nv.primitiveField()*Cw);

    forAll(N, celli)
    {
        if (N[celli] > 1e-3)
        {
            Gtot[celli] += k1_*mag(U[celli])/Lobs[celli];
        }
    }

    return tGtot;
}


Foam::tmp<Foam::volScalarField> Foam::XiGModels::basicSubGrid::Db() const
{
    const objectRegistry& db = Su_.db();
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");
    const volScalarField& rho = db.lookupObject<volScalarField>("rho");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const volScalarField& Lobs = db.lookupObject<volScalarField>("Lobs");

    return XiGModel_->Db()
        + rho*Su_*(Xi - 1.0)*mgb*(0.5*Lobs)*Lobs/(mgb*Lobs + 1.0);
}


bool Foam::XiGModels::basicSubGrid::read(const dictionary& XiGProperties)
{
    XiGModel::read(XiGProperties);
    k1_ = XiGModelCoeffs_.lookup<scalar>("k1");

    return true;
}


// ************************************************************************* //

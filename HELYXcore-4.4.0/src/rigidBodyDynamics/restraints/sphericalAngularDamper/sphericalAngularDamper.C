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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "restraints/sphericalAngularDamper/sphericalAngularDamper.H"
#include "rigidBodyModel/rigidBodyModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(sphericalAngularDamper, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        sphericalAngularDamper,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::sphericalAngularDamper::sphericalAngularDamper
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model)
{
    sphericalAngularDamper::read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::sphericalAngularDamper::~sphericalAngularDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::sphericalAngularDamper::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx
) const
{
    vector moment = -coeff_*model_.v(model_.master(bodyID_)).w();

    if (model_.debug)
    {
        Info<< " moment " << moment << endl;
    }

    // Accumulate the force for the restrained body
    fx[bodyIndex_] += spatialVector(moment, Zero);
}


bool Foam::RBD::restraints::sphericalAngularDamper::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeff_ = coeffs_.lookup<scalar>("coeff");

    return true;
}


void Foam::RBD::restraints::sphericalAngularDamper::write(Ostream& os) const
{
    restraint::write(os);

    os.writeEntry("coeff", coeff_);
}


// ************************************************************************* //

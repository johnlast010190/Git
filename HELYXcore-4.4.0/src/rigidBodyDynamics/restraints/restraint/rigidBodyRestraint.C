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

\*---------------------------------------------------------------------------*/

#include "restraints/restraint/rigidBodyRestraint.H"
#include "rigidBodyModel/rigidBodyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(restraint, 0);
    defineRunTimeSelectionTable(restraint, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraint::restraint
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    name_(name),
    bodyID_(model.bodyID(dict.lookup("body"))),
    bodyIndex_(model.master(bodyID_)),
    coeffs_(dict),
    model_(model)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraint::~restraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::RBD::restraint::coeffDict() const
{
    return coeffs_;
}


bool Foam::RBD::restraint::read(const dictionary& dict)
{
    coeffs_ = dict;
    return true;
}


void Foam::RBD::restraint::write(Ostream& os) const
{
    os.writeEntry("type", type());
    os.writeEntry("body", model_.name(bodyID_));
}


// ************************************************************************* //

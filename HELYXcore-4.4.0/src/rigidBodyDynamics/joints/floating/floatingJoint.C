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

#include "joints/floating/floatingJoint.H"
#include "rigidBodyModel/rigidBodyModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "joints/Rs/Rs.H"
#include "joints/Rzyx/Rzyx.H"
#include "joints/Pxyz/Pxyz.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(floating, 0);

    addToRunTimeSelectionTable
    (
        joint,
        floating,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBD::joints::composite>
Foam::RBD::joints::floating::sixDoF()
{
    PtrList<joint> cj(2);
    cj.set(0, new joints::Pxyz());

    // The quaternion-based spherical joint could be used
    // but then w must be set appropriately
    //cj.set(1, new joints::Rs());

    // Alternatively the Euler-angle joint can be used
    cj.set(1, new joints::Rzyx());

    return autoPtr<composite>(new composite(cj));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::floating::floating()
:
    composite(sixDoF())
{}


Foam::RBD::joints::floating::floating(const dictionary& dict)
:
    composite(sixDoF())
{}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::floating::clone() const
{
    return autoPtr<joint>(new floating(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::floating::~floating()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::floating::write(Ostream& os) const
{
    joint::write(os);
}


// ************************************************************************* //

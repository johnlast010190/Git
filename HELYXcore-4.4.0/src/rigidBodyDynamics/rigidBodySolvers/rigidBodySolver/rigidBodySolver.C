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

#include "rigidBodySolvers/rigidBodySolver/rigidBodySolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

    defineTypeNameAndDebug(rigidBodySolver, 0);
    defineRunTimeSelectionTable(rigidBodySolver, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolver::rigidBodySolver(rigidBodyMotion& body)
:
    model_(body)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolver::~rigidBodySolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBD::rigidBodySolver::correctQuaternionJoints()
{
    if (model_.unitQuaternions())
    {
        forAll(model_.joints(), i)
        {
            const label qi = model_.joints()[i].qIndex();

            if (model_.joints()[i].unitQuaternion())
            {
                // Calculate the change in the unit quaternion
                vector dv((q().block<vector>(qi) - q0().block<vector>(qi)));
                scalar magDv = mag(dv);

                if (magDv > SMALL)
                {
                    // Calculate the unit quaternion corresponding to the change
                    quaternion dQuat(dv/magDv, cos(magDv), true);

                    // Transform the previous time unit quaternion
                    quaternion quat
                    (
                        normalize(model_.joints()[i].unitQuaternion(q0())*dQuat)
                    );

                    // Update the joint unit quaternion
                    model_.joints()[i].unitQuaternion(quat, q());
                }
            }
        }
    }
}


// ************************************************************************* //

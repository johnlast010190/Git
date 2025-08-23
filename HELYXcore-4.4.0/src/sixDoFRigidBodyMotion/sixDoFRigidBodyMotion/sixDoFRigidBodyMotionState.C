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

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyMotion/sixDoFRigidBodyMotionState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState()
:
    centreOfRotation_(Zero),
    Q_(I),
    v_(Zero),
    a_(Zero),
    pi_(Zero),
    tau_(Zero)
{}


Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState
(
    const dictionary& dict
)
:
    centreOfRotation_
    (
        dict.lookupOrDefault
        (
            "centreOfRotation",
            dict.lookupOrDefault("centreOfMass", vector::zero)
        )
    ),
    Q_(dict.lookupOrDefault("orientation", tensor::I)),
    v_(dict.lookupOrDefault("velocity", vector::zero)),
    a_(dict.lookupOrDefault("acceleration", vector::zero)),
    pi_(dict.lookupOrDefault("angularMomentum", vector::zero)),
    tau_(dict.lookupOrDefault("torque", vector::zero))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionState::~sixDoFRigidBodyMotionState()
{}


// ************************************************************************* //

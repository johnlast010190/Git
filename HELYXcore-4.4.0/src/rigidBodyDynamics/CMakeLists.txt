#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
[-----------------------------------------------------------------------------]
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
    (c) 2019 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "rigidBodyDynamics")

add_helyx_library(${TARGET_NAME} SHARED
./bodies/rigidBody/rigidBody.C
./bodies/masslessBody/masslessBody.C
./bodies/jointBody/jointBody.C
./bodies/compositeBody/compositeBody.C
./bodies/subBody/subBody.C
./bodies/sphere/sphere.C
./bodies/cuboid/cuboid.C
./joints/joint/joint.C
./joints/null/nullJoint.C
./joints/composite/compositeJoint.C
./joints/floating/floatingJoint.C
./joints/Rx/Rx.C
./joints/Ry/Ry.C
./joints/Rz/Rz.C
./joints/Ra/Ra.C
./joints/Rs/Rs.C
./joints/Rzyx/Rzyx.C
./joints/Rxyz/Rxyz.C
./joints/Ryxz/Ryxz.C
./joints/Px/Px.C
./joints/Py/Py.C
./joints/Pz/Pz.C
./joints/Pa/Pa.C
./joints/Pxyz/Pxyz.C
./restraints/restraint/rigidBodyRestraint.C
./restraints/restraint/rigidBodyRestraintNew.C
./restraints/linearSpring/linearSpring.C
./restraints/linearDamper/linearDamper.C
./restraints/linearAxialAngularSpring/linearAxialAngularSpring.C
./restraints/sphericalAngularDamper/sphericalAngularDamper.C
./rigidBodyModel/rigidBodyModel.C
./rigidBodyModel/forwardDynamics.C
./rigidBodyModelState/rigidBodyModelState.C
./rigidBodyModelState/rigidBodyModelStateIO.C
./rigidBodyMotion/rigidBodyMotion.C
./rigidBodyMotion/rigidBodyMotionIO.C
./rigidBodySolvers/rigidBodySolver/rigidBodySolver.C
./rigidBodySolvers/rigidBodySolver/newRigidBodySolver.C
./rigidBodySolvers/symplectic/symplectic.C
./rigidBodySolvers/Newmark/Newmark.C
./rigidBodySolvers/CrankNicolson/CrankNicolson.C
)

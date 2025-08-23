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

set(TARGET_NAME "sixDoFRigidBodyMotion")

add_helyx_library(${TARGET_NAME} SHARED
./sixDoFRigidBodyMotion/sixDoFRigidBodyMotion.C
./sixDoFRigidBodyMotion/sixDoFRigidBodyMotionIO.C
./sixDoFRigidBodyMotion/sixDoFRigidBodyMotionState.C
./sixDoFRigidBodyMotion/sixDoFRigidBodyMotionStateIO.C
./sixDoFRigidBodyMotion/restraints/sixDoFRigidBodyMotionRestraint/sixDoFRigidBodyMotionRestraint.C
./sixDoFRigidBodyMotion/restraints/sixDoFRigidBodyMotionRestraint/sixDoFRigidBodyMotionRestraintNew.C
./sixDoFRigidBodyMotion/restraints/linearAxialAngularSpring/linearAxialAngularSpring.C
./sixDoFRigidBodyMotion/restraints/linearSpring/linearSpring.C
./sixDoFRigidBodyMotion/restraints/sphericalAngularSpring/sphericalAngularSpring.C
./sixDoFRigidBodyMotion/restraints/tabulatedAxialAngularSpring/tabulatedAxialAngularSpring.C
./sixDoFRigidBodyMotion/restraints/linearDamper/linearDamper.C
./sixDoFRigidBodyMotion/restraints/sphericalAngularDamper/sphericalAngularDamper.C
./sixDoFRigidBodyMotion/constraints/sixDoFRigidBodyMotionConstraint/sixDoFRigidBodyMotionConstraint.C
./sixDoFRigidBodyMotion/constraints/sixDoFRigidBodyMotionConstraint/sixDoFRigidBodyMotionConstraintNew.C
./sixDoFRigidBodyMotion/constraints/axis/sixDoFRigidBodyMotionAxisConstraint.C
./sixDoFRigidBodyMotion/constraints/line/sixDoFRigidBodyMotionLineConstraint.C
./sixDoFRigidBodyMotion/constraints/orientation/sixDoFRigidBodyMotionOrientationConstraint.C
./sixDoFRigidBodyMotion/constraints/plane/sixDoFRigidBodyMotionPlaneConstraint.C
./sixDoFRigidBodyMotion/constraints/point/sixDoFRigidBodyMotionPointConstraint.C
./pointPatchFields/derived/sixDoFRigidBodyDisplacement/sixDoFRigidBodyDisplacementPointPatchVectorField.C
./pointPatchFields/derived/uncoupledSixDoFRigidBodyDisplacement/uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField.C
./sixDoFRigidBodyMotionSolver/sixDoFRigidBodyMotionSolver.C
./sixDoFRigidBodyGIBMotionSolver/sixDoFRigidBodyGIBMotionSolver.C
./sixDoFSolvers/sixDoFSolver/sixDoFSolver.C
./sixDoFSolvers/sixDoFSolver/newSixDoFSolver.C
./sixDoFSolvers/symplectic/symplectic.C
./sixDoFSolvers/CrankNicolson/CrankNicolson.C
./sixDoFSolvers/Newmark/Newmark.C
)


helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
    PRIVATE surfMesh
    PRIVATE forces
    PRIVATE dynamicMesh
    PRIVATE finiteVolume
)

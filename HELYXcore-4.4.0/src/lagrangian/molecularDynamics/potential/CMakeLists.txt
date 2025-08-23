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

set(TARGET_NAME "potential")

add_helyx_library(${TARGET_NAME} SHARED
./potential/potential.C
./pairPotential/pairPotentialList/pairPotentialList.C
./pairPotential/basic/pairPotential.C
./pairPotential/basic/pairPotentialIO.C
./pairPotential/basic/pairPotentialNew.C
./pairPotential/derived/lennardJones/lennardJones.C
./pairPotential/derived/maitlandSmith/maitlandSmith.C
./pairPotential/derived/azizChen/azizChen.C
./pairPotential/derived/exponentialRepulsion/exponentialRepulsion.C
./pairPotential/derived/coulomb/coulomb.C
./pairPotential/derived/dampedCoulomb/dampedCoulomb.C
./pairPotential/derived/noInteraction/noInteraction.C
./energyScalingFunction/basic/energyScalingFunction.C
./energyScalingFunction/basic/energyScalingFunctionNew.C
./energyScalingFunction/derived/shifted/shifted.C
./energyScalingFunction/derived/shiftedForce/shiftedForce.C
./energyScalingFunction/derived/noScaling/noScaling.C
./energyScalingFunction/derived/sigmoid/sigmoid.C
./energyScalingFunction/derived/doubleSigmoid/doubleSigmoid.C
./tetherPotential/tetherPotentialList/tetherPotentialList.C
./tetherPotential/basic/tetherPotential.C
./tetherPotential/basic/tetherPotentialNew.C
./tetherPotential/derived/harmonicSpring/harmonicSpring.C
./tetherPotential/derived/restrainedHarmonicSpring/restrainedHarmonicSpring.C
./tetherPotential/derived/pitchForkRing/pitchForkRing.C
./electrostaticPotential/electrostaticPotential.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/lagrangian/molecularDynamics
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
)

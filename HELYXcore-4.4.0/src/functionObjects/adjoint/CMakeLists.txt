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

set(TARGET_NAME "adjointFunctionObjects")

add_helyx_library(${TARGET_NAME} SHARED
    ./patchBoundaryDist/patchBoundaryDist.C
    ./objectiveFunctionObject/objectiveFunctionObject.C
    ./powerLoss/powerLossObjectiveFunctionObject.C
    ./volumePowerLoss/volumePowerLossObjectiveFunctionObject.C
    ./flowRate/flowRateObjectiveFunctionObject.C
    ./uniformity/uniformityObjectiveFunctionObject.C
    ./uniformityIndex/uniformityIndexFunctionObject.C
    ./massflowSplit/massflowSplitObjectiveFunctionObject.C
    ./swirl/swirlObjectiveFunctionObject.C
    ./force/forceObjectiveFunctionObject.C
    ./rearLift/rearLiftObjectiveFunctionObject.C
    ./volumeFlowrateUniformity/volumeFlowrateUniformityObjectiveFunctionObject.C
    ./volumeUniformity/volumeUniformityObjectiveFunctionObject.C
    ./composite/compositeObjectiveFunctionObject.C
    ./shearStress/shearStressObjectiveFunctionObject.C
    ./moment/momentObjectiveFunctionObject.C
    ./pumpEfficiency/pumpEfficiencyObjectiveFunctionObject.C
    ./neutral/neutralObjectiveFunctionObject.C
    ./forceOverForce/forceOverForceObjectiveFunctionObject.C
    ./concentration/concentrationObjectiveFunctionObject.C
    ./torque/torqueObjectiveFunctionObject.C
    ./turbulentNoise/turbulentNoiseObjectiveFunctionObject.C
    ./volumeVelocityLimit/volumeVelocityLimitObjectiveFunctionObject.C
    ./temperature/temperatureObjectiveFunctionObject.C
    ./temperatureVolume/temperatureVolumeObjectiveFunctionObject.C
    ./potential/potentialObjectiveFunctionObject.C
    ./reversePotential/reversePotentialObjectiveFunctionObject.C
    ./velocityComponents/velocityComponentsObjectiveFunctionObject.C
    )

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fileFormats-includes
        ${HELYX_PROJECT_DIR}/src/transportModels
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fluidThermophysicalModels
    PRIVATE incompressibleTransportModels
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE specie
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE fvOptions
    PRIVATE eulerianMultiphaseSystems
    PRIVATE multiphaseMaterials
    PRIVATE materialModels
    PRIVATE solidThermo
    )

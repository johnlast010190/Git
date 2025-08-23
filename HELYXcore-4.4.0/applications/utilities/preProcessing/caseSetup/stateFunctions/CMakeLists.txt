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

set(TARGET_NAME "stateFunctions")

add_helyx_library(${TARGET_NAME} SHARED
./modificationSwitches/modificationSwitches.C
./stateFunction/stateFunction.C
./stateComponents/region/regionState.C
./stateComponents/turbulenceModel/turbulenceModelState.C
./stateComponents/singlePhaseFluid/singlePhaseFluidState.C
./stateComponents/solid/solidState.C
./stateComponents/multiphaseFluid/multiphaseFluidState.C
./stateComponents/multiphaseVOF/multiphaseVOFState.C
./stateComponents/multiphaseEuler/multiphaseEulerState.C
./stateComponents/multiphaseDriftFlux/multiphaseDriftFluxState.C
./stateComponents/highMachFluid/highMachFluidState.C
./stateComponents/multicomponentFluid/multicomponentFluidState.C
./stateComponents/multiphaseMulticomponent/multiphaseMulticomponentState.C
./stateComponents/VOF/VOFState.C
./stateComponents/EulerEuler/EulerEulerState.C
./stateComponents/porousSolid/porousSolidState.C
./stateFunctions/global/globalState.C
./stateFunctions/null/nullState.C
./stateFunctions/simpleFoam/simpleFoamState.C
./stateFunctions/simplecFoam/simplecFoamState.C
./stateFunctions/helyxSegregated/helyxSegregatedState.C
./stateFunctions/tHelyxSegregatedRAS/tHelyxSegregatedRASState.C
./stateFunctions/tHelyxSegregatedLES/tHelyxSegregatedLESState.C
./stateFunctions/tHelyxSegregatedVOFRAS/tHelyxSegregatedVOFRASState.C
./stateFunctions/tHelyxSegregatedVOFLES/tHelyxSegregatedVOFLESState.C
./stateFunctions/tHelyxSegregatedEulerEulerRAS/tHelyxSegregatedEulerEulerRASState.C
./stateFunctions/tHelyxSegregatedEulerEulerLES/tHelyxSegregatedEulerEulerLESState.C
./stateFunctions/helyxCoupled/helyxCoupledState.C
./stateFunctions/rhoHelyxCoupled/rhoHelyxCoupledState.C
./stateFunctions/trhoHelyxCoupledRAS/trhoHelyxCoupledRASState.C
./stateFunctions/trhoHelyxCoupledLES/trhoHelyxCoupledLESState.C
./stateFunctions/tHelyxCoupledRAS/tHelyxCoupledRASState.C
./stateFunctions/tHelyxCoupledLES/tHelyxCoupledLESState.C
./stateFunctions/rhoSimpleFoam/rhoSimpleFoamState.C
./stateFunctions/rhoSimplecFoam/rhoSimplecFoamState.C
./stateFunctions/pimpleFoamRAS/pimpleFoamRASState.C
./stateFunctions/pimpleFoamLES/pimpleFoamLESState.C
./stateFunctions/rhoPimpleFoamRAS/rhoPimpleFoamRASState.C
./stateFunctions/rhoPimpleFoamLES/rhoPimpleFoamLESState.C
./stateFunctions/rhoCentralFoamRAS/rhoCentralFoamRASState.C
./stateFunctions/rhoCentralFoamLES/rhoCentralFoamLESState.C
./stateFunctions/hisa/hisaState.C
./stateFunctions/hisaRAS/hisaRASState.C
./stateFunctions/hisaLES/hisaLESState.C
./stateFunctions/reactingFoamRAS/reactingFoamRASState.C
./stateFunctions/rhoReactingBuoyantFoamLES/rhoReactingBuoyantFoamLESState.C
./stateFunctions/buoyantMulticomponent/buoyantMulticomponentState.C
./stateFunctions/buoyantMulticomponentRAS/buoyantMulticomponentRASState.C
./stateFunctions/buoyantMulticomponentLES/buoyantMulticomponentLESState.C
./stateFunctions/porous/porousState.C
./stateFunctions/tPorous/tPorousState.C
./stateFunctions/chtSolid/chtSolidState.C
./stateFunctions/chtTransientSolid/chtTransientSolidState.C
./stateFunctions/interFoamRAS/interFoamRASState.C
./stateFunctions/interFoamLES/interFoamLESState.C
./stateFunctions/compressibleInterFoamRAS/compressibleInterFoamRASState.C
./stateFunctions/multiphaseEulerFoamRAS/multiphaseEulerFoamRASState.C
./stateFunctions/multiphaseEulerFoamLES/multiphaseEulerFoamLESState.C
./stateFunctions/multiphaseDriftFluxFoamRAS/multiphaseDriftFluxFoamRASState.C
./stateFunctions/adjointMultipoint/adjointMultipointState.C
./stateIndex/stateIndex.C
./fieldInitialisation/fieldInit/fieldInit.C
./fieldInitialisation/fieldInit/newFieldInit.C
./fieldInitialisation/boundaryValue/boundaryValue.C
./fieldInitialisation/buoyant/buoyant.C
./fieldInitialisation/cellSetMethod/cellSetMethod.C
./fieldInitialisation/fixedValue/fixedValue.C
./fieldInitialisation/Laplace/Laplace.C
./fieldInitialisation/meanSurfaceIntensity/meanSurfaceIntensity.C
./fieldInitialisation/none/none.C
./fieldInitialisation/potentialFlow/potentialFlow.C
./fieldInitialisation/Prandtl/Prandtl.C
./fieldInitialisation/relativeHumidity/relativeHumidity.C
./fieldInitialisation/turbulentIL/turbulentIL.C
./fieldInitialisation/buoyantPressure/buoyantPressure.C
./fieldInitialisation/volumeHeight/volumeHeight.C
./fieldInitialisation/equationOfState/equationOfState.C
)

helyx_link_libraries(${TARGET_NAME}
PRIVATE meshTools
PRIVATE surfMesh
PRIVATE fileFormats
PRIVATE sampling
PRIVATE surfMesh
PRIVATE finiteVolume
PRIVATE genericPatchFields
PRIVATE fluidThermophysicalModels
PRIVATE materialModels
PRIVATE multiphaseMaterials
PRIVATE solidThermo
PRIVATE OpenFOAM
PRIVATE lagrangian
PRIVATE regionModels
PRIVATE fvOptions
PRIVATE eulerianMultiphaseSystems
PRIVATE materialModels
)

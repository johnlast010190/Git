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
    (c) 2019-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "compressibleTurbulenceModels")

add_helyx_library(${TARGET_NAME} SHARED
./compressibleTurbulenceModel.C
./turbulentFluidThermoModels/turbulentFluidThermoModels.C
./turbulentFluidThermoModels/derivedFvPatchFields/turbulentConvectiveTemperature/turbulentConvectiveTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.C
./turbulentFluidThermoModels/derivedFvPatchFields/turbulentTemperatureCoupledBaffleMixed/turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/thermalBaffle1D/thermalBaffle1DFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/totalFlowRateAdvectiveDiffusive/totalFlowRateAdvectiveDiffusiveFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/turbulentTemperatureRadCoupledMixed/turbulentTemperatureRadCoupledMixedFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/externalWallHeatFluxTemperature/externalWallHeatFluxTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallHeatTransfer/wallHeatTransferFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/convectiveHeatTransfer/convectiveHeatTransferFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/fixedIncidentRadiation/fixedIncidentRadiationFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/lumpedMassWallTemperature/lumpedMassWallTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/solidNodeHeatFluxTemperature/solidNodeHeatFluxTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/lumpedServerOutletTemperature/lumpedServerOutletTemperatureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/outletMappedUniformInletHeatAddition/outletMappedUniformInletHeatAdditionFvPatchField.C
./turbulentFluidThermoModels/derivedFvPatchFields/outletMachNumberPressure/outletMachNumberPressureFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatWallFunction/alphatWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatPhaseChangeWallFunction/alphatPhaseChangeWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtWallFunction/DtWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtJayatillekeWallFunction/DtJayatillekeWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtKaderWallFunction/DtKaderWallFunctionFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/groupOfHumansHeatFlux/groupOfHumansHeatFluxFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/groupOfHumansHumidityFlux/groupOfHumansHumidityFluxFvPatchScalarField.C
./turbulentFluidThermoModels/derivedFvPatchFields/groupOfHumansCO2/groupOfHumansCO2FvPatchScalarField.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        SLGThermo-includes
    PUBLIC
        fluidThermophysicalModels-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE radiationModels
    PRIVATE materialModels
    PRIVATE multiphaseMaterials
    PRIVATE fluidThermophysicalModels
    PRIVATE radiationModels
    PRIVATE solidThermo
    PRIVATE specie
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE knn
    PRIVATE helyxMap
    PUBLIC turbulenceModels
)

helyx_link_libraries(${TARGET_NAME}-includes
    INTERFACE turbulenceModels-includes
)

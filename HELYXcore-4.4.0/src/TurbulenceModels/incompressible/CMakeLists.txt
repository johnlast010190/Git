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

set(TARGET_NAME "incompressibleTurbulenceModels")

add_helyx_library(${TARGET_NAME} SHARED
./incompressibleTurbulenceModel.C
./turbulentTransportModels/turbulentTransportModels.C
./turbulentTransportModels/RAS/qZeta/qZeta.C
./turbulentTransportModels/RAS/LamBremhorstKE/LamBremhorstKE.C
./turbulentTransportModels/RAS/LienLeschziner/LienLeschziner.C
./turbulentTransportModels/RAS/ShihQuadraticKE/ShihQuadraticKE.C
./turbulentTransportModels/RAS/LienCubicKE/LienCubicKE.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatWallFunction/alphatWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtJayatillekeWallFunction/DtJayatillekeWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtKaderWallFunction/DtKaderWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtWallFunction/DtWallFunctionFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/turbulentConvectiveTemperature/turbulentConvectiveTemperatureFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/environmentalTemperature/environmentalTemperatureFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/groupOfHumansHeatFlux/groupOfHumansHeatFluxFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/groupOfHumansHumidityFlux/groupOfHumansHumidityFluxFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/groupOfHumansCO2/groupOfHumansCO2FvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/solidNodeHeatFluxTemperature/solidNodeHeatFluxTemperatureFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/lumpedMassWallTemperature/lumpedMassWallTemperatureFvPatchScalarField.C
./turbulentTransportModels/derivedFvPatchFields/lumpedServerOutletTemperature/lumpedServerOutletTemperatureFvPatchScalarField.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fluidThermophysicalModels-includes
        specie-includes
        solidThermo-includes
)

if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    helyx_link_libraries(${TARGET_NAME}
    PRIVATE incompressibleTransportModels
    PUBLIC turbulenceModels
    PRIVATE radiationModels
    PUBLIC finiteVolume
    PRIVATE ODE
    PRIVATE meshTools
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE solidThermo
    PRIVATE knn
    PRIVATE helyxMap
    )

elseif(${HELYX_SYSTEM_NAME} STREQUAL "POSIX")
    helyx_link_libraries(${TARGET_NAME}
    PRIVATE incompressibleTransportModels
    PRIVATE turbulenceModels
    PRIVATE radiationModels
    PRIVATE finiteVolume
    PRIVATE ODE
    PRIVATE meshTools
    PRIVATE knn
    PRIVATE helyxMap
    )
endif ()

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

set(TARGET_NAME "phaseCompressibleTurbulenceModels")

add_helyx_library(${TARGET_NAME} SHARED
./phaseCompressibleTurbulenceModels.C
./phasePressureModel/phasePressureModel.C
./kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.C
./kineticTheoryModels/viscosityModel/viscosityModel/viscosityModel.C
./kineticTheoryModels/viscosityModel/viscosityModel/newViscosityModel.C
./kineticTheoryModels/viscosityModel/Gidaspow/GidaspowViscosity.C
./kineticTheoryModels/viscosityModel/Syamlal/SyamlalViscosity.C
./kineticTheoryModels/viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
./kineticTheoryModels/viscosityModel/none/noneViscosity.C
./kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.C
./kineticTheoryModels/conductivityModel/conductivityModel/newConductivityModel.C
./kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.C
./kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.C
./kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C
./kineticTheoryModels/radialModel/radialModel/radialModel.C
./kineticTheoryModels/radialModel/radialModel/newRadialModel.C
./kineticTheoryModels/radialModel/CarnahanStarling/CarnahanStarlingRadial.C
./kineticTheoryModels/radialModel/LunSavage/LunSavageRadial.C
./kineticTheoryModels/radialModel/SinclairJackson/SinclairJacksonRadial.C
./kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.C
./kineticTheoryModels/granularPressureModel/granularPressureModel/newGranularPressureModel.C
./kineticTheoryModels/granularPressureModel/Lun/LunPressure.C
./kineticTheoryModels/granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C
./kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.C
./kineticTheoryModels/frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
./kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
./kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
./kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C
./kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaFvPatchScalarField.C
./kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/transportModel
        phaseCompressibleTurbulenceModel-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE incompressibleTransportModels
    PRIVATE compressibleTwoPhaseSystem
    PRIVATE finiteVolume
    PRIVATE fvOptions
    PRIVATE meshTools
)


set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.C
    ${CMAKE_CURRENT_LIST_DIR}/kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.C
    ${CMAKE_CURRENT_LIST_DIR}/kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

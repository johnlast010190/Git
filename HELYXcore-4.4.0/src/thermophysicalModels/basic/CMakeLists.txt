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

set(TARGET_NAME "fluidThermophysicalModels")

add_helyx_library(${TARGET_NAME} SHARED
./basicThermo/basicThermo.C
./fluidThermo/fluidThermo.C
./rhoThermo/rhoThermo.C
./rhoThermo/rhoThermos.C
./addStaticHead/staticHead.C
./derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.C
./derivedFvPatchFields/fixedEnergyZone/fixedEnergyZoneFvPatchScalarField.C
./derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.C
./derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.C
./derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.C
./derivedFvPatchFields/energyJump/energyJump/energyJumpFvPatchScalarField.C
./derivedFvPatchFields/energyJump/energyJumpAMI/energyJumpAMIFvPatchScalarField.C
./derivedFvPatchFields/coupledTotalPressure/coupledTotalPressureFvPatchScalarField.C
./derivedFvPatchFields/totalPressure/totalPressureFvPatchScalarField.C
./derivedFvPatchFields/timeVaryingMappedTotalPressure/timeVaryingMappedTotalPressureFvPatchScalarField.C
./derivedFvPatchFields/fanPressure/fanPressureFvPatchScalarField.C
./derivedFvPatchFields/fanPressureAMI/fanPressureAMIFvPatchScalarField.C
./derivedFvPatchFields/targetFanPressurePoint/targetFanPressurePointOutletPressureFvPatchScalarField.C
./derivedFvPatchFields/outletInletPressure/outletInletPressureFvPatchScalarField.C
./derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.C
./derivedFvPatchFields/pressureVelocity/pressureVelocityFvPatchVectorField.C
./derivedFvPatchFields/rotatingTotalPressure/rotatingTotalPressureFvPatchScalarField.C
./derivedFvPatchFields/fixedprghPressure/fixedprghPressureFvPatchScalarField.C
./planarAverage/planarAverage.C

# ddt models compatible with material models
./finiteVolume/ddtSchemes/smartBoundedDdtScheme/smartBoundedDdtSchemes.C

./porosityModels/DarcyForchheimer/DarcyForchheimer.C
./porosityModels/DarcyForchheimer/temperatureDependentDarcyForchheimer/temperatureDependentDarcyForchheimer.C

# Reacting thermo part
./mixtures/basicSpecieMixture/basicSpecieMixture.C
./mixtures/basicCombustionMixture/basicCombustionMixture.C
./fluidMulticomponentThermo/fluidMulticomponentThermo.C
./psiuMulticomponentThermo/psiuMulticomponentThermo.C
./psiuMulticomponentThermo/psiuMulticomponentThermos.C
./rhoMulticomponentThermo/rhoMulticomponentThermo.C
./rhoMulticomponentThermo/rhoMulticomponentThermos.C
./derivedFvPatchFields/fixedUnburntEnthalpy/fixedUnburntEnthalpyFvPatchScalarField.C
./derivedFvPatchFields/gradientUnburntEnthalpy/gradientUnburntEnthalpyFvPatchScalarField.C
./derivedFvPatchFields/mixedUnburntEnthalpy/mixedUnburntEnthalpyFvPatchScalarField.C
./functionObjects/moleFractions/moleFractions.C
./functionObjects/massFractions/massFractions.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        meshTools-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE specie
    PUBLIC thermophysicalMixtures
    PRIVATE materialModels
    PRIVATE thermophysicalProperties
    PRIVATE finiteVolume
    PRIVATE finiteArea
    PRIVATE meshTools
)

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

set(TARGET_NAME "radiationModels")

add_helyx_library(${TARGET_NAME} SHARED
./radiationModels/radiationModel/radiationModel.C
./radiationModels/radiationModel/radiationModelNew.C
./radiationModels/noRadiation/noRadiation.C
./radiationModels/P1/P1.C
./radiationModels/fvDOM/fvDOM/fvDOM.C
./radiationModels/fvDOM/radiativeIntensityRay/radiativeIntensityRay.C
./radiationModels/fvDOM/blackBodyEmission/blackBodyEmission.C
./radiationModels/fvDOM/absorptionCoeffs/absorptionCoeffs.C
./radiationModels/viewFactor/viewFactor.C
./radiationModels/opaqueSolid/opaqueSolid.C
./radiationModels/solarLoad/solarLoad.C
./radiationModels/solarLoad/faceShading/faceShading.C
./radiationModels/domSolar/domSolar.C
./radiationModels/domSolar/radiantField.C
./scatterModels/scatterModel/scatterModel.C
./scatterModels/scatterModel/scatterModelNew.C
./scatterModels/noScatter/noScatter.C
./scatterModels/constantScatter/constantScatter.C
./absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.C
./absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModelNew.C
./absorptionEmissionModels/multiBandAbsorptionEmission/multiBandAbsorptionEmission.C
./absorptionEmissionModels/multiBandZoneAbsorptionEmission/multiBandZoneAbsorptionEmission.C
./absorptionEmissionModels/noAbsorptionEmission/noAbsorptionEmission.C
./absorptionEmissionModels/constantAbsorptionEmission/constantAbsorptionEmission.C
./absorptionEmissionModels/binary/binary.C
./absorptionEmissionModels/greyMean/greyMean.C
./absorptionEmissionModels/wideBand/wideBand.C
./absorptionEmissionModels/greyMeanSolid/greyMeanSolid.C
./absorptionEmissionModels/multiBandSolid/multiBandSolid.C
./absorptionEmissionModels/TaylorFoster/TaylorFoster.C
./boundaryRadiationProperties/boundaryRadiationProperties.C
./boundaryRadiationProperties/boundaryRadiationPropertiesPatch.C
./sootModels/sootModel/sootModel.C
./sootModels/sootModel/sootModelNew.C
./sootModels/noSoot/noSoot.C
./transmissivityModels/transmissivityModel/transmissivityModel.C
./transmissivityModels/transmissivityModel/transmissivityModelNew.C
./transmissivityModels/noTransmissivity/noTransmissivity.C
./transmissivityModels/constantTransmissivity/constantTransmissivity.C
./transmissivityModels/multiBandSolidTransmissivity/multiBandSolidTransmissivity.C
./solarCalculator/solarCalculator.C
./derivedFvPatchFields/MarshakRadiation/MarshakRadiationFvPatchScalarField.C
./derivedFvPatchFields/MarshakRadiationFixedTemperature/MarshakRadiationFixedTemperatureFvPatchScalarField.C
./derivedFvPatchFields/greyDiffusiveRadiation/greyDiffusiveRadiationMixedFvPatchScalarField.C
./derivedFvPatchFields/wideBandDiffusiveRadiation/wideBandDiffusiveRadiationMixedFvPatchScalarField.C
./derivedFvPatchFields/greyDiffusiveViewFactor/greyDiffusiveViewFactorFixedValueFvPatchScalarField.C
)


helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        surfMesh-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE solidThermo
    PRIVATE SLGThermo
    PRIVATE thermophysicalProperties
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE distributed
    PRIVATE materialModels
    PRIVATE multiphaseMaterials
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/transmissivityModel/noTransmissivity/noTransmissivity.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )

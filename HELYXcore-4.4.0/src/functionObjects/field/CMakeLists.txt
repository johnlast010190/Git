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

set(TARGET_NAME "fieldFunctionObjects")

add_helyx_library(${TARGET_NAME} SHARED
./columnAverage/columnAverage.C
./fieldAverage/fieldAverage.C
./fieldAverage/fieldAverageItem/fieldAverageItem.C
./fieldAverage/fieldAverageItem/fieldAverageItemIO.C
./reverseFieldAverage/reverseFieldAverage.C
./reverseFieldAverage/reverseFieldAverageItem/reverseFieldAverageItem.C
./reverseFieldAverage/reverseFieldAverageItem/reverseFieldAverageItemIO.C
./fieldCoordinateSystemTransform/fieldCoordinateSystemTransform.C
./fieldMinMax/fieldMinMax.C
./fieldValues/fieldValue/fieldValue.C
./fieldValues/fieldValue/fieldValueNew.C
./fieldValues/fieldValueDelta/fieldValueDelta.C
./fieldValues/volFieldValue/volFieldValue.C
./fieldValues/surfaceFieldValue/surfaceFieldValue.C
./fieldCalcs/fieldCalc/fieldCalc.C
./fieldCalcs/fieldCalc/fieldCalcNew.C
./fieldCalcs/surfaceFieldCalc/surfaceFieldCalc.C
./nearWallFields/nearWallFields.C
./nearWallFields/findCellParticle.C
./nearWallFields/findCellParticleCloud.C
./processorField/processorField.C
./readFields/readFields.C
./setFlow/setFlow.C
./streamLine/streamLine.C
./streamLine/streamLineBase.C
./streamLine/streamLineParticle.C
./streamLine/streamLineParticleCloud.C
./wallBoundedStreamLine/wallBoundedStreamLine.C
./wallBoundedStreamLine/wallBoundedStreamLineParticle.C
./wallBoundedStreamLine/wallBoundedStreamLineParticleCloud.C
./wallBoundedStreamLine/wallBoundedParticle.C
./surfaceInterpolate/surfaceInterpolate.C
./surfaceReconstruct/surfaceReconstruct.C
./regionSizeDistribution/regionSizeDistribution.C
./histogram/histogram.C
./fieldExpression/fieldExpression.C
./components/components.C
./randomise/randomise.C
./Curle/Curle.C
./div/div.C
./grad/grad.C
./mag/mag.C
./magSqr/magSqr.C
./vorticity/vorticity.C
./enstrophy/enstrophy.C
./Q/Q.C
./Qdot/Qdot.C
./Lambda2/Lambda2.C
./flowType/flowType.C
./CourantNo/CourantNo.C
./aniCourantNo/aniCourantNo.C
./PecletNo/PecletNo.C
./blendingFactor/blendingFactor.C
./pressure/pressure.C
./MachNo/MachNo.C
./visibilityIndex/visibilityIndex.C
./scale/scale.C
./fieldsExpression/fieldsExpression.C
./add/add.C
./subtract/subtract.C
./turbulenceFields/turbulenceFields.C
./yPlus/yPlus.C
./wallShearStress/wallShearStress.C
./wallHeatFlux/wallHeatFlux.C
./writePatchDist/writePatchDist.C
./writeCellCentres/writeCellCentres.C
./writeCellVolumes/writeCellVolumes.C
./XiReactionRate/XiReactionRate.C
./streamFunction/streamFunction.C
./valueAverage/valueAverage.C
./statisticalQuantities/statisticalQuantities.C
./fluxSummary/fluxSummary.C
./mapFields/mapFields.C
./reactionSensitivityAnalysis/reactionsSensitivityAnalysisObjects.C
./DESModelRegions/DESModelRegions.C
./externalCoupled/externalCoupled.C
./externalCoupled/externalCoupledMixed/externalCoupledMixedFvPatchFields.C
./externalCoupled/externalCoupledTemperatureMixed/externalCoupledTemperatureMixedFvPatchScalarField.C
./extractEulerianParticles/extractEulerianParticles/extractEulerianParticles.C
./extractEulerianParticles/eulerianParticle/eulerianParticle.C
./particleDistribution/particleDistribution.C
./flux/flux.C
./ddt2/ddt2.C
./zeroGradient/zeroGradient.C
./fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.C
./fieldBlendingFactor/blendingFunctions/heaviside/heaviside.C
./fieldBlendingFactor/blendingFunctions/cubicRamp/cubicRamp.C
./fieldBlendingFactor/blendingFunctions/linearRamp/linearRamp.C
./fieldBlendingFactor/blendingFunctions/tanhRamp/tanhRamp.C
./fieldBlendingFactor/blendingSources/blendingSource/blendingSource.C
./fieldBlendingFactor/blendingSources/cellSizeRatioBlendingSource/cellSizeRatioBlendingSource.C
./fieldBlendingFactor/blendingSources/courantNumberBlendingSource/courantNumberBlendingSource.C
./fieldBlendingFactor/blendingSources/pecletBlendingSource/pecletBlendingSource.C
./fieldBlendingFactor/blendingSources/phaseBlendingSource/phaseBlendingSource.C
./fieldBlendingFactor/blendingSources/faceSetBlendingSource/faceSetBlendingSource.C
./fieldBlendingFactor/blendingSources/zonalDESBlendingSource/zonalDESBlendingSource.C
./fieldBlendingFactor/blendingSources/nearWallBlendingSource/nearWallBlendingSource.C
./fieldBlendingFactor/blendingSources/hybridBlendingSource/hybridBlendingSource.C
./fieldBlendingFactor/blendingSources/cellGradientBlendingSource/cellGradientBlendingSource.C
./fieldBlendingFactor/blendingSources/fieldBlendingSource/fieldBlendingSource.C
./fieldBlendingFactor/blendingSources/residualBlendingSource/residualBlendingSource.C
./fieldBlendingFactor/fieldBlendingItem/fieldBlendingItem.C
./fieldBlendingFactor/fieldBlendingFactor/fieldBlendingFactor.C
./fieldProcess/fieldProcess.C
./fieldProcess/fieldProcessFunctions.C
./fieldWrite/fieldWrite.C
./proudmanAcoustics/proudmanAcoustics.C
./fractionalEffectiveDose/fractionalEffectiveDose.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fileFormats-includes
        ${HELYX_PROJECT_DIR}/src/transportModels
        twoPhaseMixture-includes
        specie-includes
        materialModels-includes
        twoPhaseMixture-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE lagrangian
    PRIVATE distributionModels
    PRIVATE sampling
    PRIVATE surfMesh
    PRIVATE materialModels
    PRIVATE fluidThermophysicalModels
    PRIVATE incompressibleTransportModels
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE combustionModels
    PRIVATE chemistryModel
    PRIVATE solidThermo
    PRIVATE pairPatchAgglomeration
    PRIVATE fvOptions
    PRIVATE eulerianMultiphaseSystems
    PRIVATE multiphaseMaterials
    PRIVATE materialModels
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/wallBoundedStreamLine/wallBoundedStreamLine.C
    ${CMAKE_CURRENT_LIST_DIR}/wallBoundedStreamLine/wallBoundedStreamLineParticleCloud.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

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
    (c) 2019-2021 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "waves2FoamProcessing")

add_helyx_library(${TARGET_NAME} SHARED
./preProcessing/setWaveProperties/setWaveProperties.C
./preProcessing/setWaveProperties/current/potentialCurrentProperties/potentialCurrentProperties.C
./preProcessing/setWaveProperties/solitary/solitaryFirstProperties/solitaryFirstProperties.C
./preProcessing/setWaveProperties/regular/stokesFirstProperties/stokesFirstProperties.C
./preProcessing/setWaveProperties/regular/stokesFirstStandingProperties/stokesFirstStandingProperties.C
./preProcessing/setWaveProperties/regular/stokesSecondProperties/stokesSecondProperties.C
./preProcessing/setWaveProperties/regular/stokesSecondModulationProperties/stokesSecondModulationProperties.C
./preProcessing/setWaveProperties/regular/stokesFifthProperties/stokesFifthProperties.C
./preProcessing/setWaveProperties/irregular/bichromaticFirstProperties/bichromaticFirstProperties.C
./preProcessing/setWaveProperties/irregular/bichromaticSecondProperties/bichromaticSecondProperties.C
./preProcessing/setWaveProperties/irregular/waveSpectra/waveSpectra.C
./preProcessing/setWaveProperties/irregular/waveSpectra/frequencyAxis/frequencyAxis.C
./preProcessing/setWaveProperties/irregular/waveSpectra/frequencyAxis/equidistantFrequencyAxis/equidistantFrequencyAxis.C
./preProcessing/setWaveProperties/irregular/waveSpectra/frequencyAxis/cosineStretchedFrequencyAxis/cosineStretchedFrequencyAxis.C
./preProcessing/setWaveProperties/irregular/waveSpectra/spectra/JONSWAP/JONSWAP.C
./preProcessing/setWaveProperties/irregular/waveSpectra/spectra/PiersonMoskowitz/PiersonMoskowitz.C
./preProcessing/setWaveProperties/irregular/waveSpectra/phases/phases.C
./preProcessing/setWaveProperties/irregular/waveSpectra/phases/randomPhase/randomPhase.C
./preProcessing/setWaveProperties/irregular/waveSpectra/phases/focusingPhase/focusingPhase.C
./preProcessing/setWaveProperties/irregular/irregularProperties/irregularProperties.C
./preProcessing/setWaveProperties/mixedType/combinedWavesProperties/combinedWavesProperties.C
./preProcessing/setWaveProperties/externalWaveForcing/externalSourceProperties/externalSourceProperties.C
./preProcessing/setWaveProperties/externalWaveForcing/oceanWave3DProperties/oceanWave3DProperties.C
./preProcessing/probes/waveGauges/waveGauges.C
./preProcessing/probes/probeGauges/probeGauges.C
./preProcessing/probes/pointDistributions/pointDistributions.C
./preProcessing/probes/pointDistributions/userDefinedDistribution/userDefinedDistribution.C
./preProcessing/probes/pointDistributions/lineDistribution/lineDistribution.C
./preProcessing/probes/pointDistributions/circularDistribution/circularDistribution.C
./preProcessing/probes/pointDistributions/quadrilateralDistribution/quadrilateralDistribution.C
./postProcessing/postProcessingWaves/postProcessingWaves.C
./postProcessing/postProcessingWaves/writeRawData/rawSurfaceElevation/rawSurfaceElevation.C
./postProcessing/postProcessingWaves/writeRawData/rawVelocityProbes/rawVelocityProbes.C
./postProcessing/postProcessingWaves/writeRawData/rawForcesAndMoments/rawForcesAndMoments.C
./postProcessing/postProcessingWaves/writeRawData/rawAlphaProbes/rawAlphaProbes.C
./postProcessing/postProcessingWaves/writeRawData/rawOvertopping/rawOvertopping.C
./postProcessing/postProcessingWaves/write2Ascii/write2Ascii.C
./postProcessing/postProcessingWaves/writeIndexLocation/writeIndexLocation.C
./postProcessing/postProcessingWaves/cleaningUtil/removeData/removeData.C
./postProcessing/postProcessingWaves/cleaningUtil/removeSubDirs/removeSubDirs.C
./postProcessing/postProcessingWaves/writeInterpolatedData/interpolateSurfaceElevation/interpolateSurfaceElevation.C
./postProcessing/postProcessingWaves/writeInterpolatedData/interpolateVelocityProbes/interpolateVelocityProbes.C
./postProcessing/postProcessingWaves/writeInterpolatedData/interpolateForcesAndMoments/interpolateForcesAndMoments.C
./postProcessing/postProcessingWaves/writeInterpolatedData/interpolateAlphaProbes/interpolateAlphaProbes.C
./postProcessing/postProcessingWaves/writeInterpolatedData/interpolateOvertopping/interpolateOvertopping.C
./postProcessing/postProcessingWaves/spectralAnalysis/spectralMethods/complexExp.C
./postProcessing/postProcessingWaves/spectralAnalysis/spectralMethods/spectralMethodsLeastSquaresBased.C
./postProcessing/postProcessingWaves/spectralAnalysis/leastSquareBasedMethods/powerSpectraLS/powerSpectraLS.C
./postProcessing/postProcessingWaves/spectralAnalysis/leastSquareBasedMethods/reflectionAnalysis2DLS/reflectionAnalysis2DLS.C
./postProcessing/postProcessingWaves/miscOperations/ensembleAverage/ensembleAverage.C
./postProcessing/postProcessingWaves/miscOperations/trapz/trapz.C
./postProcessing/postProcessingWaves/miscOperations/cumTrapz/cumTrapz.C
./postProcessing/postProcessingWaves/miscOperations/zeroCrossing/zeroCrossing.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        lagrangian-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE sampling
    PRIVATE fileFormats
    PRIVATE waves2Foam
)

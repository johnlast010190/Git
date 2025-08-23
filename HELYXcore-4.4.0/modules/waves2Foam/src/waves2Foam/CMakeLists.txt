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

set(TARGET_NAME "waves2Foam")

add_helyx_library(${TARGET_NAME} SHARED
    ./convexPolyhedral/convexPolyhedral.C
    ./convexPolyhedral/localFace.C
    ./convexPolyhedral/localCell.C
    ./convexPolyhedral/localCellNeg.C
    ./waveTheories/waveTheory/waveTheory.C
    ./waveTheories/waveTheory/newWaveTheory.C
    ./waveTheories/current/potentialCurrent/potentialCurrent.C
    ./waveTheories/solitary/solitaryFirst/solitaryFirst.C
    ./waveTheories/regular/stokesFirst/stokesFirst.C
    ./waveTheories/regular/stokesFirstStanding/stokesFirstStanding.C
    ./waveTheories/regular/stokesSecond/stokesSecond.C
    ./waveTheories/regular/stokesSecondModulation/stokesSecondModulation.C
    ./waveTheories/regular/stokesFifth/stokesFifth.C
    ./waveTheories/regular/streamFunction/streamFunction.C
    ./waveTheories/irregular/irregular/irregular.C
    ./waveTheories/irregular/bichromaticFirst/bichromaticFirst.C
    ./waveTheories/irregular/bichromaticSecond/bichromaticSecond.C
    ./waveTheories/mixedType/combinedWaves/combinedWaves.C
    ./fvPatchFields/waveAlpha/waveAlphaFvPatchScalarField.C
    ./fvPatchFields/wavePressure/wavePressureFvPatchScalarField.C
    ./fvPatchFields/waveVelocity/waveVelocityFvPatchVectorField.C
    ./fvPatchFields/wavePressureC/wavePressureCFvPatchScalarField.C
    ./fvPatchFields/waveVelocityC/waveVelocityCFvPatchVectorField.C
    ./relaxationZone/relaxationShape/relaxationShape.C
    ./relaxationZone/relaxationShape/newRelaxationShape.C
    ./relaxationZone/relaxationShape/frozen/relaxationShapeFrozen.C
    ./relaxationZone/relaxationShape/rectangular/relaxationShapeRectangular.C
    ./relaxationZone/relaxationShape/cylindrical/relaxationShapeCylindrical.C
    ./relaxationZone/relaxationShape/cylindricalMoving/relaxationShapeCylindricalMoving.C
    ./relaxationZone/relaxationShape/semiCylindrical/relaxationShapeSemiCylindrical.C
    ./relaxationZone/relaxationShape/line/relaxationShapeLine.C
    ./relaxationZone/relaxationWeight/relaxationWeight.C
    ./relaxationZone/relaxationWeight/relaxationWeightCourantCorrection.C
    ./relaxationZone/relaxationWeight/exponential/relaxationWeightExponential.C
    ./relaxationZone/relaxationWeight/freePolynomial/relaxationWeightFreePolynomial.C
    ./relaxationZone/relaxationWeight/thirdOrderPolynomial/relaxationWeightThirdOrderPolynomial.C
    ./relaxationZone/relaxationScheme/relaxationScheme.C
    ./relaxationZone/relaxationScheme/newRelaxationScheme.C
    ./relaxationZone/relaxationScheme/empty/relaxationSchemeEmpty.C
    ./relaxationZone/relaxationScheme/spatial/relaxationSchemeSpatial.C
    ./relaxationZone/relaxationScheme/spatialInterpolation/relaxationSchemeSpatialInterpolation.C
    ./relaxationZone/relaxationScheme/spatialImplicit/relaxationSchemeSpatialImplicit.C
    ./relaxationZone/numericalBeach/numericalBeach.C
    ./relaxationZone/numericalBeach/newNumericalBeach.C
    ./relaxationZone/numericalBeach/empty/numericalBeachEmpty.C
    ./relaxationZone/relaxationZone.C
    ./setWaveField/setWaveField.C
    ./include/crossVersionCompatibility.C
    ./solverObjects/waves/waves.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE sampling
    PRIVATE fileFormats
    PRIVATE fluidThermophysicalModels
    PRIVATE multiphaseMaterials
    PRIVATE solidThermo
    PRIVATE materialModels
    PRIVATE helyxSolve
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        lagrangian-includes
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/relaxationZone/relaxationScheme/spatial/relaxationSchemeSpatial.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

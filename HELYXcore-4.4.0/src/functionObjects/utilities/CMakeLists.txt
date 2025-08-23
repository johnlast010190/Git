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
    (c) 2019-2020 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "utilityFunctionObjects")

add_helyx_library(${TARGET_NAME} SHARED
./abort/abort.C
./codedFunctionObject/codedFunctionObject.C
./combineMultipointFields/combineMultipointFields.C
./ensightWrite/ensightWrite.C
./ensightWrite/ensightWriteUpdate.C
./vtkWrite/vtkWrite.C
./vtkWrite/vtkWriteUpdate.C
./removeRegisteredObject/removeRegisteredObject.C
./residuals/residuals.C
./runTimeControl/runTimeControl.C
./runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.C
./runTimeControl/runTimeCondition/runTimeCondition/runTimeConditionNew.C
./runTimeControl/runTimeCondition/equationMaxIterCondition/equationMaxIterCondition.C
./runTimeControl/runTimeCondition/equationInitialResidualCondition/equationInitialResidualCondition.C
./runTimeControl/runTimeCondition/minMaxCondition/minMaxCondition.C
./runTimeControl/runTimeCondition/averageCondition/averageCondition.C
./runTimeControl/runTimeCondition/minTimeStepCondition/minTimeStepCondition.C
./runTimeControl/runTimeCondition/minIterCondition/minIterCondition.C
./runTimeControl/runTimeCondition/confidenceIntervalCondition/confidenceIntervalCondition.C
./runTimeControl/runTimeCondition/reverseAnalysisStatisticalCondition/reverseAnalysisStatisticalCondition.C
./setTimeStep/setTimeStepFunctionObject.C
./systemCall/systemCall.C
./timeActivatedFileUpdate/timeActivatedFileUpdate.C
./writeDictionary/writeDictionary.C
./writeObjects/writeObjects.C
./writeLastTime/writeLastTime.C
./thermoCoupleProbes/thermoCoupleProbes.C
./thermalComfort/thermalComfort.C
./thermalComfort/comfortFunctions.C
./acousticAnalogy/FfowcsWilliamsHawkings.C
./acousticAnalogy/DavidsonCurle.C
./acousticAnalogy/soundObserver.C
./acousticAnalogy/kfwhModel.C
./acousticAnalogy/solidObject.C
./acousticAnalogy/observer/observer.C
./acousticAnalogy/observer/patchObserver/patchObserver.C
./acousticAnalogy/observer/setObserver/setObserver.C
./acousticAnalogy/observer/vibroacousticObserver/vibroacousticObserver.C
./inSituFourier/inSituFourier.C
./twoPointCorrelation/twoPointCorrelation.C
./heatBalance/heatBalance.C
./writeSensor/writeSensor.C
./pedestrianComfort/pedestrianComfort.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/../thirdParty/knn/src
        fileFormats-includes
        ${HELYX_PROJECT_DIR}/src/transportModels
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE conversion
    PRIVATE sampling
    PRIVATE fluidThermophysicalModels
    PRIVATE radiationModels
    PRIVATE incompressibleTransportModels
    PRIVATE ODE
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE dynamicMesh
    PRIVATE knn
    PRIVATE helyxMap
    PRIVATE helyxSolve
)

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

set(TARGET_NAME "report")

add_helyx_library(${TARGET_NAME} SHARED
./convergenceTermination/convergenceTermination.C
./subVolume/subVolume.C
./fieldInstance/fieldInstance.C
./fieldInstance/fieldInstanceIO.C
./volumeReport/volumeReport.C
./flowCharacteristics/flowCharacteristics.C
./ccnvReport/ccnvReport.C
./courant/courant.C
./fluidHeight/fluidHeight.C
./patchWrite/patchWrite.C
./printDictionary/printDictionary.C
./surfaceReport/surfaceReport.C
./surfaceAreaReport/surfaceAreaReport.C
./massInOutReport/massInOutReport.C
./waveGauge/waveGauge.C
./waveRunup/waveRunup.C
./centerOfGravity/centerOfGravity.C
./wallContactTime/wallContactTime.C
./catchRatio/catchRatio.C
./highSpeedFlowOutletMonitor/highSpeedFlowOutletMonitor.C
./framesMonitor/framesMonitor.C
)



helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE fileFormats
    PRIVATE dynamicMesh
    PRIVATE dynamicFvMesh
    PRIVATE sampling
    PRIVATE fluidThermophysicalModels
    PRIVATE materialModels
)

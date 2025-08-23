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

set(TARGET_NAME "helyxMap")

add_helyx_library(${TARGET_NAME} SHARED
./containers/Array2d/Array2d.C
./containers/Array3d/Array3d.C
./containers/dataItem/dataItem.C
./knn/knnInterpolation.C
./meshes/smartBoundBox/smartBoundBox.C
./sample/fieldSample/fieldSample.C
./sample/gridField/gridField.C
./sample/helyxSample/helyxSample.C
./sample/helyxSamplePar/helyxSamplePar.C
./sample/npGrid/npGrid.C
./sample/trainGeoModel/trainGeoModel.C
./surfaceMap/surfaceMap.C
./tools/mapScalar/mapScalar.C
./tools/readFields/readFields.C
./helyxMap.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        surfMesh-includes
)

target_include_directories(${TARGET_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${HELYX_PROJECT_DIR}/thirdParty/knn/src>
        $<INSTALL_INTERFACE:thirdParty/knn/src>
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE knn
)

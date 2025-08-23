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

set(TARGET_NAME "fileFormats")

helyx_flex_target(STLReaderASCII
    ./stl/STLReaderASCII.L
    ${CMAKE_CURRENT_BINARY_DIR}/STLReaderASCII.C
    COMPILE_FLAGS "-+ -f"
    )

add_helyx_library(${TARGET_NAME} SHARED
./ensight/file/ensightCase.C
./ensight/file/ensightCaseOptions.C
./ensight/file/ensightFile.C
./ensight/file/ensightGeoFile.C
./ensight/part/ensightCells.C
./ensight/part/ensightFaces.C
./ensight/read/ensightReadFile.C
./ensight/type/ensightPTraits.C
./nas/NASCore.C
./csv/CSVCore.C
./fire/FIRECore.C
./starcd/STARCDCore.C
./stl/STLCore.C
./stl/STLReader.C
${FLEX_STLReaderASCII_OUTPUTS}
./vtk/file/foamVtkFileWriter.C
./vtk/file/foamVtkSeriesWriter.C
./vtk/file/foamVtmWriter.C
./vtk/core/foamVtkCore.C
./vtk/core/foamVtkPTraits.C
./vtk/format/foamVtkFormatter.C
./vtk/format/foamVtkAsciiFormatter.C
./vtk/format/foamVtkBase64Formatter.C
./vtk/format/foamVtkAppendBase64Formatter.C
./vtk/format/foamVtkAppendRawFormatter.C
./vtk/format/foamVtkBase64Layer.C
./vtk/format/foamVtkLegacyAsciiFormatter.C
./vtk/format/foamVtkLegacyRawFormatter.C
./vtk/output/foamVtkOutput.C
./vtk/output/foamVtkOutputOptions.C
./vtk/part/foamVtkMeshMaps.C
./vtk/part/foamVtuCells.C
./vtk/part/foamVtuSizing.C
./vtk/read/vtkUnstructuredReader.C
./vtk/write/foamVtkSurfaceWriter.C
./coordSet/coordSet.C
./f4st/read/f4stReader.C
./sampledSetWriters/writers.C
./sampledSetWriters/ensight/ensightSetWriterRunTime.C
./sampledSetWriters/gnuplot/gnuplotSetWriterRunTime.C
./sampledSetWriters/jplot/jplotSetWriterRunTime.C
./sampledSetWriters/raw/rawSetWriterRunTime.C
./sampledSetWriters/vtk/vtkSetWriterRunTime.C
./sampledSetWriters/xmgrace/xmgraceSetWriterRunTime.C
./sampledSetWriters/csv/csvSetWriterRunTime.C
./sampledSetWriters/tsv/tsvSetWriterRunTime.C
./sampledSetWriters/f4st/f4stSetWriterRunTime.C
)


# Let CMake deals with Flex headers unless FLEX_ARCH_PATH is set.
# Make sure the right Flex headers are included in case FLEX_ARCH_PATH is set
# and a different Flex version exist in /usr/*
#
# N.B.: Looks like there is a bug in CMake find Flex module when using tool-chains/windows:
# Flex executable is found and the right includes were used but
# FLEX_INCLUDE_DIR is reported as not found (FLEX_INCLUDE_DIR-NOTFOUND)
# Therefore, only include the header if FLEX_ARCH_PATH is set for linux
if(NOT "" STREQUAL "${FLEX_ARCH_PATH}" AND NOT "${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    helyx_additional_includes(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_FLEX_INC}
    )
endif()

if("${OPENCASCADE_FOUND}")
    target_sources(${TARGET_NAME}
        PRIVATE ./cad/CADReader.C
    )

    helyx_link_libraries(${TARGET_NAME}
        PUBLIC ${THIRDPARTY_OPENCASCADE}
    )

    set_source_files_properties(
        ${CMAKE_CURRENT_LIST_DIR}/cad/CADReader.C
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )
endif()

if("${ALEMBIC_FOUND}")
    target_sources(${TARGET_NAME}
        PRIVATE ./alembic/AlembicReader.C
    )

    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_ALEMBIC}
    )

    set_source_files_properties(
        ${CMAKE_CURRENT_LIST_DIR}/alembic/AlembicReader.C
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )
endif()

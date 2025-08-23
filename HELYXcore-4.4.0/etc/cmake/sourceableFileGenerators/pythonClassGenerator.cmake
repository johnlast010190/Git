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

Description
    Generates a sourceable Python file

[----------------------------------------------------------------------------]]

get_filename_component(REAL_TP_DIR "${HELYX_THIRDPARTY_DIR}" REALPATH)
file(RELATIVE_PATH RELATIVE_TP_DIR ${HELYX_PROJECT_DIR} ${HELYX_THIRDPARTY_DIR})
string(REPLACE
    "ThirdParty-${HELYX_THIRDPARTY_VERSION}"
    'ThirdParty-' + self.HELYX_THIRDPARTY_VERSION
    TP_DIR_STRING
    ${RELATIVE_TP_DIR}
    )
set(HELYX_THIRDPARTY_DIR_STRING "os.path.join(self.HELYX_PROJECT_DIR, 'ThirdParty-' + self.HELYX_THIRDPARTY_VERSION)")


if("${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/bin")
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "os.path.join(self.HELYX_PROJECT_DIR, 'platforms', self.HELYX_OPTIONS, 'bin')")
else()
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}\"")
endif()

if("${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib")
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "os.path.join(self.HELYX_PROJECT_DIR, 'platforms', self.HELYX_OPTIONS, 'lib')")
else()
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "\"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}\"")
endif()


if(NOT "" STREQUAL "${MPI_ARCH_PATH}")
    get_filename_component(REAL_MPI_DIR "${MPI_ARCH_PATH}" REALPATH)
    string(REPLACE
        "${REAL_TP_DIR}/" ""
        MPI_DIR_STRING "${REAL_MPI_DIR}"
        )
    if(NOT "${MPI_DIR_STRING}" STREQUAL "${REAL_MPI_DIR}")
        set(MPI_DIR_STRING "os.path.join(self.HELYX_THIRDPARTY_DIR, '${MPI_DIR_STRING}')")
    endif()
else()
    set(REAL_MPI_DIR "")
    set(MPI_DIR_STRING "\"\"")
endif()

# Always populate extraSourceableFiles dir to make packing reproducible
# As a transition, save both '${HELYX_OPTIONS}.py' and '${HELYX_OPTIONS}_${HELYX_MPI_NAME}.py'
set(PYTHON_FILE ${CMAKE_BINARY_DIR}/extraSourceableFiles/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.py)

configure_file("${CMAKE_CURRENT_LIST_DIR}/pythonClassTemplate.py" "${PYTHON_FILE}"
    @ONLY
    NEWLINE_STYLE UNIX
)

if ("python" IN_LIST RUNTIME_SHELL)
    file(COPY "${PYTHON_FILE}"
        DESTINATION "${HELYX_PROJECT_DIR}/platforms/"
    )
endif()

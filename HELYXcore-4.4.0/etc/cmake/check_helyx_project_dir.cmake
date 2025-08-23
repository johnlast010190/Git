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
    (c) 2019-2022 Engys Ltd.

Description
    A cmake script used to check if HELYX_PROJECT_DIR is properly set,
    otherwise, error.

[----------------------------------------------------------------------------]]

# ============================================================================ #
# -------------------------- Check HELYX_PROJECT_DIR ------------------------- #
# ============================================================================ #

# N.B.  This need to be a file to be included in userSettings.cmake
#        Proobably we will never reach any error from this file

# HELYX_PROJECT_DIR is defined on userSettings.cmake,
# so it make sense to double check is HELYX_PROJECT_DIR is a valid path

# HELYX_PROJECT_DIR is usually defined automatically, which shouldn`t cause any issue.
# However, it can be manually defined elsewhere and therefore
# recieve any string, which could be potentially wrong

# TODO: Need to check for symlinks?


# First, ensure that only one path is stored in HELYX_PROJECT_DIR
list(LENGTH HELYX_PROJECT_DIR helyx_project_dir_size)
if(NOT "${HELYX_PROJECT_DIR}" STREQUAL "" AND helyx_project_dir_size GREATER 1 )

    string(CONCAT s
        "HELYX_PROJECT_DIR should contain only one string. "
        "However, multiple strings were found:\n"
        "\"${HELYX_PROJECT_DIR}\"\n"
        )
    message(FATAL_ERROR "${s}")

elseif(NOT "${HELYX_PROJECT_DIR}" STREQUAL "" AND IS_DIRECTORY "${HELYX_PROJECT_DIR}")

    # Now, check if HELYX_PROJECT_DIR and CMAKE_CURRENT_LIST_DIR are consistent
    get_filename_component(cmake_current_project "${CMAKE_CURRENT_LIST_DIR}/../../" REALPATH)
    get_filename_component(expected_helyx_project "${HELYX_PROJECT_DIR}" REALPATH)

    if(NOT "${cmake_current_project}" MATCHES "${expected_helyx_project}")
        string(CONCAT s
            "HELYX_PROJECT_DIR is inconsistent with the current project directory. "
            "Please, review the HELYX setup in the userSettings file.\n"
            "CURRENT PROJECT DIR: \"${cmake_current_project}\"\n"
            "HELYX_PROJECT_DIR:   \"${expected_helyx_project}\"\n"
            )
        message(FATAL_ERROR "${s}")
    endif()

    # Finally, check if it is a valid HELYX project based on CMake
    # I think we can use CMakeLists.txt as a north to define if
    # it is a valid HELYX prooject
    if(EXISTS "${HELYX_PROJECT_DIR}/CMakeLists.txt")

        file(STRINGS "${HELYX_PROJECT_DIR}/CMakeLists.txt" valid_helyx_project REGEX "^project\\(HELYXcore\\)")

        if(NOT valid_helyx_project)
            string(CONCAT s
                "HELYX_PROJECT_DIR evaluated as the following path "
                "\"${HELYX_PROJECT_DIR}\"\n"
                "Is not a HELYX project based on CMake. "
                "Please, review the HELYX setup in the userSettings file.\n"
            )
            message(FATAL_ERROR "${s}")
        endif()

    else()
        string(CONCAT s
            "HELYX_PROJECT_DIR evaluated as the following string, which doesn`t "
            "seem to be a valid HELYX project:\n"
            "\"${HELYX_PROJECT_DIR}\"\n"
            "Please, review the HELYX setup in the userSettings file.\n"
        )
        message(FATAL_ERROR "${s}")
    endif()

else()
    string(CONCAT s
        "HELYX_PROJECT_DIR evaluated as the following string, which is not a "
        "valid path:\n"
        "\"${HELYX_PROJECT_DIR}\"\n"
    )
    message(FATAL_ERROR "${s}")
endif()

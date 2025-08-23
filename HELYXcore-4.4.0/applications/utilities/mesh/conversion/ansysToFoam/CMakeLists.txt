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

set(TARGET_NAME "ansysToFoam-exe")

helyx_flex_target(ansysToFoam
    ./ansysToFoam.L
    ${CMAKE_CURRENT_BINARY_DIR}/ansysToFoam.C
    COMPILE_FLAGS "-+ -f"
    )

add_helyx_executable(${TARGET_NAME}
${FLEX_ansysToFoam_OUTPUTS}
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

helyx_link_libraries(${TARGET_NAME}
    PRIVATE OpenFOAM
)

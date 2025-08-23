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
    (c) 2021 Engys Ltd

Description
    Script used to check Ninja version when -clean* flags are used in emake
    in this case, ninja-1.10 give spurious errors related to lnInclude dirs

[----------------------------------------------------------------------------]]

execute_process(
    COMMAND ninja --version
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    OUTPUT_VARIABLE NINJA_VERSION
)

# check if NINJA_VERSION matches the X.X.X pattern
if(NOT NINJA_VERSION MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    unset(NINJA_VERSION)
endif()

if(NINJA_VERSION)
    STRING(CONCAT s
        "\nNinja may leads to spurious errors "
        "when using emake with \"-clean\" flags.\n "
        "These Ninja errors will not affect HELYX build."
        )
    message(WARNING "${s}")
endif()
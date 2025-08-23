#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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

if ("${HELYX_PROJECT_DIR}" STREQUAL "")
    message(FATAL_ERROR "Environment has not been set; cannot continue.")
endif()


add_subdirectory(calcNonUniformOffsetsForMapped)

add_subdirectory(fieldReport)

add_subdirectory(funkyDoCalc)

add_subdirectory(funkySetBoundaryField)

add_subdirectory(funkySetFields)

add_subdirectory(funkySetLagrangianField)

add_subdirectory(funkyWarpMesh)

add_subdirectory(replayTransientBC)

add_subdirectory(writeBoundarySubfields)

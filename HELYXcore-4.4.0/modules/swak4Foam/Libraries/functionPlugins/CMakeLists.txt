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


add_subdirectory(lagrangianCloudAdaptors)

add_subdirectory(swakChemistryModelFunctionPlugin)

add_subdirectory(swakDynamicMeshFunctionPlugin)

add_subdirectory(swakFvcSchemesFunctionPlugin)

add_subdirectory(swakLagrangianCloudSourcesFunctionPlugin)

add_subdirectory(swakLocalCalculationsFunctionPlugin)

add_subdirectory(swakMeshQualityFunctionPlugin)

add_subdirectory(swakMeshWaveFunctionPlugin)

add_subdirectory(swakMRFFunctionPlugin)

add_subdirectory(swakQuantileFunctionPlugin)

add_subdirectory(swakRadiationModelFunctionPlugin)

add_subdirectory(swakRandomFunctionPlugin)

add_subdirectory(swakShiftFieldFunctionPlugin)

add_subdirectory(swakSurfacesAndSetsFunctionPlugin)

add_subdirectory(swakThermoTurbFunctionPlugin)

add_subdirectory(swakTransportTurbFunctionPlugin)

add_subdirectory(swakVelocityFunctionPlugin)

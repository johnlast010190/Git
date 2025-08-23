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
    (c) 2019 Engys Ltd

[----------------------------------------------------------------------------]]

set(TARGET_NAME "swakStateMachine")

add_helyx_library(${TARGET_NAME} SHARED
./BaseClasses/StateMachine.C
./BaseClasses/StateMachineRepository.C
./FunctionObjects/stateMachineCreateAndUpdateFunctionObject.C
./FunctionObjects/stateMachineSetStateFunctionObject.C
./FunctionObjects/stateMachineStateFunctionObject.C
./PluginFunctions/isStateFunctionPlugin.C
./PluginFunctions/timeSinceChangeFunctionPlugin.C
./PluginFunctions/stepsSinceChangeFunctionPlugin.C
./PluginFunctions/changedToFunctionPlugin.C
./other/stateMachineFvSolutionFvSchemesFunctionObject.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        surfMesh-includes
        fileFormats-includes
        finiteArea-includes
        sampling-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE swakLagrangianParser
    PRIVATE simpleFunctionObjects
    PRIVATE meshTools
    PRIVATE specie
    PRIVATE finiteVolume
)
helyx_link_libraries_with_generated_headers(${TARGET_NAME} PRIVATE swak4FoamParsers)

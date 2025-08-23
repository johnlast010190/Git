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

set(TARGET_NAME "reactingTwoPhaseSystem")

add_helyx_library(${TARGET_NAME} SHARED
./diameterModels/IATE/IATE.C
./diameterModels/IATE/IATEsources/IATEsource/IATEsource.C
./diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.C
./diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.C
./diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.C
./diameterModels/IATE/IATEsources/phaseChange/phaseChange.C
./twoPhaseSystem.C
./newTwoPhaseSystem.C
./twoPhaseSystems.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fluidThermophysicalModels-includes
        specie-includes
        phaseCompressibleTurbulenceModel-includes
        turbulenceModels-includes
        finiteVolume-includes
        meshTools-includes
        sampling-includes
)

helyx_link_libraries(${TARGET_NAME}
PRIVATE compressibleTurbulenceModels
PRIVATE combustionModels
PRIVATE reactingPhaseSystem
)

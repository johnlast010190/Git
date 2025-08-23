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

[----------------------------------------------------------------------------]]

set(TARGET_NAME "hydroStaticFoam-exe")

add_helyx_executable(${TARGET_NAME}
    ./hydroStaticFoam.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE immiscibleIncompressibleTwoPhaseMixture
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE finiteVolume
    PRIVATE fvOptions
    PRIVATE meshTools
    PRIVATE sampling
    PRIVATE waves2Foam
    PRIVATE waves2FoamSampling
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        twoPhaseMixture-includes
        interfaceProperties-includes
        incompressibleTransportModels-includes
        ${HELYX_SRC}/transportModels  # Should be PUBLIC in immiscibleIncompressibleTwoPhaseMixture (or something)
)

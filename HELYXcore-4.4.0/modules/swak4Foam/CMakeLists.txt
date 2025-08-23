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
    (c) 2019-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

# Swak fails to compile without BISON, so BISON_REQUIRED should default to ON
# unless otherwise set (note that this cache variable doesn't use FORCE).
set(BISON_REQUIRED ON CACHE STRING "")

# If it reaches this file, is because swak4Foam_FOUND is TRUE

# BISON is mandatory for swak4Foam compilation
# so, first set some sensible variables based on possible user combinations of
# ${module}_REQUIRED and ${dependency}_REQURIED
helyx_check_module_mandatory_dependency(swak4Foam BISON)

# if swak4Foam compilation is disabled, simply return from this script
if(swak4Foam_REQUIRED STREQUAL "OFF" OR NOT swak4Foam_FOUND)
    return()
endif()

# At this stage, BISON_REQUIRED should be ON

# check whether BISON is available
helyx_add_module_dependency(swak4Foam BISON)

# Remove BISON-generated warnings from HELYX compilation
add_compile_options(
-Wno-deprecated-copy
-Wno-deprecated-declarations
-Wno-unused-but-set-variable
)

include("etc/cmake/functions.cmake")

add_subdirectory(Libraries)

add_subdirectory(Utilities)

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
    (c) 2019-2021 Engys Ltd

[----------------------------------------------------------------------------]]

message(TITLE "Module configuration")
message(CLEAN "")

set(modules
    waves2Foam;
    HELYX-marine;
    HELYX-coupled;
    HELYX-adjoint;
    runTimePostProcessing;
    swak4Foam
)

# If the user has specified additional modules, they probably want configuration
# to fail if those modules aren't found.  We could set REQUIRED variables to ON
# by default, but probably better to be consistent.
# foreach(additional_module IN LISTS HELYX_ADDITIONAL_MODULES)
#     set(${additional_module}_REQUIRED ON
#         CACHE STRING
#         "Trinary switch to toggle the optional compilation of the \"${additional_module}\" module"
#     )
# endforeach()

list(APPEND modules "${HELYX_ADDITIONAL_MODULES}")

# Note that the module variable is used to differentiate module binaries, so
# some behaviour (particularly packing) is dependent on "module" taking the
# correct value when add_subdirectory is called.  That's probably not a great
# idea.
foreach(module ${modules})

    initialise_search_variables(${module})

    # There's an argument to be made that the COMPILE_ variables controlled
    # whether HELYX compiled modules, and the _REQUIRED variables controlled
    # whether HELYX compiled *against* ThirdParty libraries.  However, it's
    # easier to remember how to turn off module compilation if it shares the
    # same pattern as ThirdParty (particularly because modules are rarely
    # turned off, whereas ThirdParty libs are all the time).
    # Ultimately, practicality beats purity.
    if(DEFINED COMPILE_${module} OR NOT "" STREQUAL "${COMPILE_${module}}")
        string(CONCAT s
        "COMPILE_${module} was set to \"${COMPILE_${module}}\".\n"
        "The \"COMPILE_<module>\" variables were replaced with variables of "
        "the form \"<module>_REQUIRED\" in HELYX 3.4, in order to be "
        "consistent with similar variables used for ThirdParty library "
        "detection.\n"
        "Please replace any references to \"COMPILE_${module}\" in your user "
        "settings file with \"${module}_REQUIRED\".\n"
        )
        message(WARNING "${s}")
    endif()

    # Now start searching
    if(NOT "${${module}_REQUIRED}" MATCHES "^ON$|^AUTO$")
        set(${module}_FOUND_MESSAGE
            "${module}_REQUIRED set to ${${module}_REQUIRED}"
            CACHE INTERNAL ""
        )
    else()
        # ${module}_REQUIRED is ON or AUTO, so search
        if(EXISTS "${HELYX_PROJECT_DIR}/modules/${module}/CMakeLists.txt")
            set(${module}_FOUND TRUE CACHE INTERNAL "")

            # Each module gets its own findModule directory
            # This goes some way to enforcing consistency for module
            # dependencies
            set(${module}_FINDMODULE_DIR
                "${HELYX_PROJECT_DIR}/modules/${module}/etc/cmake/findModules"
                CACHE STRING
                "The directory in which to look for CMake findModules (if any) for dependencies of ${module}"
            )

            # Immediately parse the module.  This gives the module a chance to
            # declare itself missing if necessary.
            add_subdirectory(${module})
        else()
            # ${module}_REQUIRED is ON or AUTO, but module directory doesn't
            # exist
            set(${module}_FOUND_MESSAGE
                "Nothing found at \"${HELYX_PROJECT_DIR}/modules/${module}\""
                CACHE INTERNAL ""
            )
            set(${module}_FOUND FALSE CACHE INTERNAL "")
        endif()
    endif()
endforeach()

# Having given modules a chance to declare themselves missing, report!
set(missing_modules "${modules}")
foreach(module ${modules})
    if(${module}_FOUND)
        list(APPEND added_modules ${module})
        list(REMOVE_ITEM missing_modules "${module}")
    endif()
endforeach()

if(NOT "${added_modules}" STREQUAL "")
    message(CLEAN "")
    message(SUBTITLE "Found modules")
    foreach(module ${added_modules})
        message(CLEAN "${module}")
        list(REMOVE_ITEM missing_modules "${module}")
    endforeach()
endif()

if(NOT "${missing_modules}" STREQUAL "")
    message(CLEAN "")
    message(SUBTITLE "Missing modules")
    foreach(module ${missing_modules})
        if("${${module}_REQUIRED}" STREQUAL "ON")
            string(CONCAT s
                "${module}\n"
                "\t${module}_REQUIRED set to ${${module}_REQUIRED}, but "
                "${module} was not found.\n"
                "\t${${module}_FOUND_MESSAGE}\n"
                )
            message(SEND_ERROR "${s}")
        else()
            message(CLEAN "${module}  (${${module}_FOUND_MESSAGE})")
        endif()
    endforeach()
endif()

message(CLEAN
"
[==============================================================================]

")

# For use elsewhere in the code
set(FOUND_MODULES "${added_modules}" CACHE INTERNAL "")

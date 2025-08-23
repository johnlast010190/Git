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
    (c) 2019-2021 Engys Ltd.

Description
    Include file for the top-level CMakeLists.txt in a project. Can be used
    in a custom project outside the source tree.

[----------------------------------------------------------------------------]]



# --------- Runtime loading and renaming for all Windows executables --------- #

if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_EXE_TARGETS)
    foreach (TARGET_NAME IN LISTS TARGETS)
        process_windows_executable(${TARGET_NAME})
    endforeach ()
endif()


# ---------------- Disabling of Unity for specific targets ------------------- #

foreach(tgt ${TARGETS_TO_EXCLUDE_FROM_UNITY})
    if(TARGET "${tgt}")
        message(STATUS "Disabling Unity build for \"${tgt}\"")
        set_target_properties("${tgt}" PROPERTIES UNITY_BUILD OFF)
    else()
        string(CONCAT s
            "Failed to turn off Unity build for target \"${tgt}\", as it is "
            "not a target!"
            )
        message(WARNING "${s}")
    endif()
endforeach()



# ------------------------- Warn if links unresolved ------------------------- #

set(unresolved_includes "${UNRESOLVED_LINKS_TO_INCLUDES}")
list(REMOVE_DUPLICATES unresolved_includes)
list(REMOVE_ITEM unresolved_includes "")
foreach(include ${unresolved_includes})
    string(CONCAT s
        "helyx_additional_includes() called with the following argument:\n"
        "\t\"${include}\"\n"
        "This looks like it should be an INTERFACE target for specifying "
        "include directories, but no target with that name exists.\n"
        "This \"target\" was included in the following file(s):"
    )
    set(files "${${include}-included-from-here}")
    list(REMOVE_ITEM files "")
    foreach(location ${files})
        string(CONCAT s
            "${s}\n"
            "\t${location}"
        )
    endforeach()
    string(CONCAT s
        "${s}\n"
        # "(If left unchecked, these would be link errors).\n"
        "The most common cause of this bug is typos and misspellings.\n"
    )
    message(SEND_ERROR "${s}")
endforeach()


# ------------------------- Report errors/warnings --------------------------- #

if(NOT "0" STREQUAL "${WARNING_COUNT}" OR NOT "0" STREQUAL "${ERROR_COUNT}")

    set(string_to_print "  Problems detected!  ")
    pad_string_with_char("${string_to_print}" "*" 80 centre)
    message(CLEAN "\n\n${string_to_print}")

    if(NOT "0" STREQUAL "${WARNING_COUNT}")
        if("1" STREQUAL "${WARNING_COUNT}")
            message(CLEAN
            "${ColourYellow}Detected ${WARNING_COUNT} warning${ColourReset}")
        else()
            message(CLEAN
            "${ColourYellow}Detected ${WARNING_COUNT} warnings${ColourReset}")
        endif()
        message(CLEAN
            "    Your configuration may not be what you're expecting."
        )
        message(CLEAN
            "    It is good practice to remove all warnings before continuing."
        )
    endif()

    if(NOT "0" STREQUAL "${ERROR_COUNT}")
        if("1" STREQUAL "${ERROR_COUNT}")
            message(CLEAN
            "${ColourBoldRed}Detected ${ERROR_COUNT} error${ColourReset}")
        else()
            message(CLEAN
            "${ColourBoldRed}Detected ${ERROR_COUNT} errors${ColourReset}")
        endif()
        message(CLEAN
            "    CMake has failed to configure successfully."
        )
        message(CLEAN
            "    You must remove all errors before continuing."
        )
    endif()
    message(CLEAN
        "\nDetails of the detected problems can be found in the above text."
        )
    message(CLEAN
        "Please read the above text carefully, change your user settings file"
        )
    message(CLEAN
        "accordingly, and refresh the cache (\"emake -r\")."
        )

    set(string_to_print "*")
    pad_string_with_char("${string_to_print}" "*" 80 right)
    message(CLEAN ${string_to_print}\n\n\n)
endif()

# Hint: you don't want this.
if (HELYX_GENERATE_LNINCLUDES)
    execute_process(COMMAND ${CMAKE_SOURCE_DIR}/bin/mkLnIncludes.sh)
endif()

set(EMAKE_BUILD FALSE
  CACHE INTERNAL
  "Force EMAKE_BUILD to FALSE. The TRUE state now only lasts until the end of the configuration, forcing it to be set to ON every time."
  FORCE # Strictly speaking unnecessary with INTERNAL variables
)

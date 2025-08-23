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
    (c) 2019-2023 Engys Ltd.

Description
    Defines custom CMake functions

[----------------------------------------------------------------------------]]


# ============================================================================ #
# ---------------------------- General utilities ----------------------------- #
# ============================================================================ #

set(WARNING_COUNT 0
    CACHE INTERNAL
    ""
    )
set(ERROR_COUNT 0
    CACHE INTERNAL
    ""
    )
macro(message MessageType)
    if("${ARGV0}" STREQUAL "CLEAN")
        _message(NOTICE "${ARGN}")
    elseif("${ARGV0}" STREQUAL "PIPE_TERMINATED")
        string(REGEX MATCHALL "${Esc}\\[[0-9;]*m"
            matches "${ARGN}")
        string(LENGTH "${matches}" matches_len)
        if(matches)
            math(EXPR pad_length "79 + (${matches_len} - 1)")
        else()
            set(pad_length 79)
        endif()
        pad_string_with_char("${ARGN}" " " "${pad_length}")
        set(string_to_print "${string_to_print}|")
        message(CLEAN "${string_to_print}")
    elseif("${ARGV0}" STREQUAL "TITLE")
        set(string_to_print " ${ARGN} ")
        pad_string_with_char("${string_to_print}" "=" 78 centre)
        set(string_to_print "[${string_to_print}]")
        message(CLEAN "${string_to_print}")
    elseif("${ARGV0}" STREQUAL "SUBTITLE")
        set(string_to_print " ${ARGN} ")
        pad_string_with_char("${string_to_print}" "-" 78 centre)
        set(string_to_print "[${string_to_print}]")
        message(CLEAN "${string_to_print}")
    else()
        if("${ARGV0}" MATCHES "^WARNING$|^SEND_ERROR$|^FATAL_ERROR$")
            set(original_location_string "Message originates at ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE}\n")
        endif()
        set(string_to_print "${original_location_string}${ARGN}")

        if("${ARGV0}" MATCHES "^FATAL_ERROR$")
            _message(${MessageType} "${ColourBoldMagenta}${string_to_print}${ColourReset}")
        elseif("${ARGV0}" MATCHES "^SEND_ERROR$")
            _message(${MessageType} "${ColourBoldRed}${string_to_print}${ColourReset}")
        elseif("${ARGV0}" MATCHES "^WARNING$|^AUTHOR_WARNING$|^CHECK_FAIL$")
            _message(${MessageType} "${ColourYellow}${string_to_print}${ColourReset}")
        elseif("${ARGV0}" MATCHES "^STATUS$|^CHECK_START$|^CHECK_PASS$")
            _message(${MessageType} "${ColourReset}${string_to_print}${ColourReset}")
        else()
            _message("${ColourBlue}${MessageType}${string_to_print}${ColourReset}")
        endif()

        set(original_location_string "")

        if("${ARGV0}" STREQUAL "WARNING")
            math(EXPR temp "${WARNING_COUNT} + 1")
            set(WARNING_COUNT ${temp}
                CACHE INTERNAL
                ""
                )
        elseif("${ARGV0}" STREQUAL "SEND_ERROR")
            math(EXPR temp "${ERROR_COUNT} + 1")
            set(ERROR_COUNT ${temp}
                CACHE INTERNAL
                ""
                )
        endif()
    endif()
endmacro()


macro(set_target_optimisation target optimisation)
    if(NOT ${HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL})
        string(CONCAT s
        "set_target_optimisation() called on target \"${target}\" in the "
        "following directory:\n\t\"${CMAKE_CURRENT_SOURCE_DIR}\",\n but "
        "HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL set to "
        "\"${HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL}\".\nSet "
        "HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL to ON or TRUE to enable "
        "set_target_optimisation().\n"
        )
        message(AUTHOR_WARNING "${s}")
    else()
        if("${optimisation}" IN_LIST CMAKE_CONFIGURATION_TYPES)
            string(TOUPPER ${optimisation} optimisation_upper)
            separate_arguments(FLAGS_LIST UNIX_COMMAND ${CMAKE_CXX_FLAGS_${optimisation_upper}})
            target_compile_options(${target}
                PRIVATE ${FLAGS_LIST}
            )
        else()
            message(SEND_ERROR
            "Optimisation level \"${optimisation}\" not recognised.  Valid optimisations are as follows:
            ${CMAKE_CONFIGURATION_TYPES}")
        endif()
    endif()
endmacro()


macro(set_optimisation optimisation)
    if(NOT ${HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL})
        string(CONCAT s
        "set_optimisation() called in the following directory:\n\t"
        "\"${CMAKE_CURRENT_SOURCE_DIR}\",\n but "
        "HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL set to "
        "\"${HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL}\".\nSet "
        "HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL to ON or TRUE to enable "
        "set_optimisation().\n"
        )
        message(AUTHOR_WARNING "${s}")
    else()
        message(STATUS
        "Setting optimisation of the following directory to \"${optimisation}\":
        ${CMAKE_CURRENT_SOURCE_DIR}"
        )

        if("${optimisation}" IN_LIST CMAKE_CONFIGURATION_TYPES)
            string(TOUPPER ${optimisation} optimisation_upper)
            add_compile_options(${CMAKE_CXX_FLAGS_${optimisation_upper}})
        else()
            message(SEND_ERROR
            "Optimisation level \"${optimisation}\" not recognised.  Valid optimisations are as follows:
            ${CMAKE_CONFIGURATION_TYPES}")
        endif()
    endif()
endmacro()


function(thirdParty_package_deprecation_warning package_name removal_version)
    if(NOT ${ARGC} EQUAL 2)
        string(CONCAT s
            "thirdParty_package_deprecation_warning() requires exactly two "
            "arguments, but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    endif()

    foreach(INDICATOR ${package_name}_ARCH_PATH ${package_name}_REQUIRED ${package_name}_PATH_SUFFIXES)
        if(DEFINED "${INDICATOR}")
            set(${INDICATOR}_STRING "${INDICATOR}, ")
            set(REMOVAL_WARNING TRUE)
        endif()
    endforeach()
    if(REMOVAL_WARNING)
        list(APPEND missing_libraries ${package_name})
        string(CONCAT ${package_name}_FOUND_MESSAGE
            "${package_name}\n"
            "\t${package_name} support was removed in HELYX version "
            "${removal_version}.\n"
            "\tPlease remove the following variables from your user settings "
            "file:\n"
            "\t\t${${package_name}_ARCH_PATH_STRING}${${package_name}_REQUIRED_STRING}${${package_name}_PATH_SUFFIXES_STRING}"
        )
        if("${${package_name}_REQUIRED}" STREQUAL "ON")
            string(CONCAT ${package_name}_FOUND_MESSAGE
                "${package_name}\n"
                "${package_name}_REQUIRED was set to ON, but ${package_name} "
                "is no longer supported and so this requirement cannot be "
                "satisfied.\n"
                "Please remove all references to \"${package_name}\" from your "
                "user settings file."
            )
        endif()
    endif()

    # As of CMake 3.11, can't have multi-line variables in the cache
    # set(${package_name}_FOUND_MESSAGE "${${package_name}_FOUND_MESSAGE}"
        # CACHE INTERNAL ""
    # )
    set(${package_name}_FOUND_MESSAGE "${${package_name}_FOUND_MESSAGE}"
        PARENT_SCOPE
    )
    set(missing_libraries "${missing_libraries}" PARENT_SCOPE)
    unset(REMOVAL_WARNING)
endfunction()


macro(add_dummy_target target_name missing_dependency)
    set(unsupported_libs FFTW)
    if(${missing_dependency} IN_LIST unsupported_libs)
        message(STATUS "Skipping ${target_name} (${missing_dependency} not supported)")
        add_custom_target(${target_name}
            COMMAND ${CMAKE_COMMAND} -E echo "The \\\"${target_name}\\\" target has been disabled, as the \\\"${missing_dependency}\\\" library is not supported when compiling with CMake."
            COMMAND ${CMAKE_COMMAND} -E echo "Please contact Engys \\(support@engys.com\\) if you believe this library should be supported."
        )
    else()
        message(STATUS "Skipping ${target_name} (${missing_dependency} not found)")
        add_custom_target(${target_name}
            COMMAND ${CMAKE_COMMAND} -E echo "The \\\"${target_name}\\\" target has been disabled, as the \\\"${missing_dependency}\\\" dependency is missing."
            COMMAND ${CMAKE_COMMAND} -E echo "To compile this target, please resolve the missing dependency, refresh the CMake cache, and try again."
        )
    endif()
    set_target_properties(${target_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
endmacro()


function(strip_exe_from_target_name TARGET_NAME EXE_NAME_VAR)
    # Remove trailing "-exe" from target name
    string(REGEX REPLACE "(.*)-exe$" "\\1" EXE_NAME ${TARGET_NAME})
    if (NOT EXE_NAME)
        set(EXE_NAME ${TARGET_NAME})
    endif ()
    set(${EXE_NAME_VAR} ${EXE_NAME} PARENT_SCOPE)
endfunction ()


# Deals with boilerplate for trinary switches:
# - Set cache variable for ${variable_name}_REQUIRED (unforced, so user-set
#   value not overwritten),
# - Checks value is valid (if invalid, resets to AUTO),
# - Sets property to ON, OFF, or AUTO,
function (initialise_search_variables base_name)
    set(${base_name}_REQUIRED AUTO
        CACHE STRING
        "Trinary switch (ON, OFF, or AUTO)"
    )
    if(NOT "${${base_name}_REQUIRED}" MATCHES "^ON$|^OFF$|^AUTO$|^$")
        message(WARNING
        "Resetting \"${base_name}_REQUIRED\" to AUTO, as it is set to unrecognised value \"${${base_name}_REQUIRED}\".  Allowed values are ON, OFF, or AUTO.")
        set(${base_name}_REQUIRED AUTO
            CACHE STRING
            "Trinary switch (ON, OFF, or AUTO)"
            FORCE
            )
    endif()
    set_property(CACHE ${base_name}_REQUIRED PROPERTY STRINGS ON OFF AUTO)

    set(${base_name}_FOUND FALSE CACHE INTERNAL "")

endfunction()


function(helyx_add_module_dependency module dependency)
    if(NOT ${ARGC} EQUAL 2)
        string(CONCAT s
            "helyx_add_module_dependency() requires exactly two arguments, "
            "but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    elseif("${module}" STREQUAL "")
        message(FATAL_ERROR "\"module\" set to \"\"")
    elseif("${dependency}" STREQUAL "")
        message(FATAL_ERROR "\"dependency\" set to \"\"")
    endif()

    if(NOT DEFINED ${module}_FINDMODULE_DIR)
        message(FATAL_ERROR "${module}_FINDMODULE_DIR was not defined")
    elseif(NOT EXISTS "${${module}_FINDMODULE_DIR}")
        string(CONCAT s
            "${module}_FINDMODULE_DIR set to the following value, which is not "
            "a directory:\n"
            "\"${${module}_FINDMODULE_DIR}\"
            "
        )
        message(FATAL_ERROR "${s}")
    endif()

    if(NOT "${${dependency}_REQUIRED}" MATCHES "ON|OFF|AUTO")
        string(CONCAT s
            "${dependency}_REQUIRED set to \"${${dependency}_REQUIRED}\", "
            "which is invalid.  ${dependency}_REQUIRED Must be one of ON, OFF, "
            "or AUTO.\n"
            "Setting ${dependency}_FOUND to FALSE.\n"
            "Setting ${module}_FOUND to FALSE.\n"
        )
        message(SEND_ERROR "${s}")
        set(${module}_FOUND FALSE CACHE INTERNAL "")
        set(${dependency}_FOUND FALSE CACHE INTERNAL "")
        set(${module}_FOUND_MESSAGE
            "${dependency}_REQUIRED set to \"${${dependency}_REQUIRED}\", which is invalid.  Must be one of ON, OFF, or AUTO"
            CACHE INTERNAL
            ""
        )
        return()
    elseif("${${dependency}_REQUIRED}" STREQUAL "OFF")
        string(CONCAT s
        "${dependency}_REQUIRED set to OFF.\n"
        "${dependency} is required to compile ${module}, so we are setting "
        "${module}_FOUND to FALSE.\n"
        "Please either set ${module}_REQUIRED to OFF (to disable ${module}), "
        "or ${dependency}_REQUIRED to \"ON\" or \"AUTO\" (to enable "
        "${dependency} detection, and therefore ${module} compilation).
        "
        )
        message(WARNING "${s}")
        # TODO:  Bug here when dependency is OFF, but module can still compile?
        set(${module}_FOUND FALSE CACHE INTERNAL "")
        set(${module}_FOUND_MESSAGE
            "Missing dependency \"${dependency}\", because ${dependency}_REQUIRED set to OFF"
            CACHE INTERNAL
            ""
        )
        return()
    endif()

    # We now know that ${module}_FINDMODULE_DIR takes a sane value, and that
    # ${module}_REQUIRED is ON or AUTO, so we search for the module.


    # Push CMAKE_MODULE_PATH (to stop to stop findModules from this HELYX module
    # infecting other parts of the code), then find dependency, then pop.
    set(pushed_module_path "${CMAKE_MODULE_PATH}")
    list(APPEND CMAKE_MODULE_PATH "${${module}_FINDMODULE_DIR}")
    find_package(${dependency})
    set(CMAKE_MODULE_PATH "${pushed_module_path}")

    # Now check whether the dependency has been found
    # Complication:  What about the case where a depndency is not required, but
    # if it's found then the module gets some extra toys?

    # dependency required to compile module:
    # - set REQUIRED to ON in userSettings.cmake
    # - set extra module REQUIRED to whatever you like (although FOUND will be
    #   false if dependency is missing)

    # Dependency just gives module extra toys:
    # - set dependency REQUIRED to AUTO (or even OFF) in user settings
    # - Set extra module REQUIRED to whatever (FOUND will be true if module
    #   found, regardless of whether dependency is found)

    # So here, just set FOUND to false if REQUIRED is ON

    if(NOT ${dependency}_FOUND)
        if(${dependency}_REQUIRED STREQUAL "ON")
            string(CONCAT s
                "${dependency} not found, but it is required to compile ${module}.\n"
                "Setting ${module}_FOUND to FALSE.\n"
                "If you do not wish to compile ${module}, please set "
                "${module}_REQUIRED to OFF.\n"
            )
            message(SEND_ERROR "${s}")
            set(${module}_FOUND FALSE CACHE INTERNAL "")
            set(${module}_FOUND_MESSAGE
                "Missing dependency \"${dependency}\""
                CACHE INTERNAL
                ""
            )
        elseif(${dependency}_REQUIRED STREQUAL "AUTO")
            string(CONCAT s
                "${dependency} not found.\n\t"
                "${dependency}_REQUIRED set to \"${${dependency}_REQUIRED}\", "
                "so ${dependency} is not required to compile ${module}.\n\t"
                "Leaving ${module}_FOUND set to \"${${module}_FOUND}\".\n\t"
                "If ${module} is compiled, some functionality may be missing.\n"
            )
            message(STATUS "${s}")
        else()
            # Shouldn't ever get here
            string(CONCAT s
                "Unknown error when finding ${dependency}.\n"
                "${dependency}_REQUIRED set to \"${${dependency}_REQUIRED}\", "
                "Should be one of ON, OFF, or AUTO.\n"
                "Setting ${module}_FOUND to FALSE.\n"
            )
            message(SEND_ERROR "${s}")
            set(${module}_FOUND FALSE CACHE INTERNAL "")
            set(${module}_FOUND_MESSAGE
                "Missing dependency \"${dependency}\""
                CACHE INTERNAL
                ""
            )
        endif()
    endif()

endfunction()

# function to be called inside modules/${module}/CMakeLists.txt
# in this way, ${module}_FOUND is TRUE if the function is called
function(helyx_check_module_mandatory_dependency module dependency)
    if(NOT ${ARGC} EQUAL 2)
        string(CONCAT s
            "helyx_check_module_mandatory_dependency() requires exactly two arguments, "
            "but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    elseif("${module}" STREQUAL "")
        message(FATAL_ERROR "\"module\" set to \"\"")
    elseif("${dependency}" STREQUAL "")
        message(FATAL_ERROR "\"dependency\" set to \"\"")
    endif()


    STRING(CONCAT message_intro
        "${module} module is found.\n"
        "${module}_REQUIRED is set to ${${module}_REQUIRED} and\ "
        "${dependency}_REQUIRED is set to ${${dependency}_REQUIRED}.\n"
        )

    if(${module}_REQUIRED STREQUAL "OFF")

        # if ${module}_REQUIRED is OFF and ${dependency}_REQUIRED is OFF:
        #   there is nothing to do

        # if ${module}_REQUIRED is OFF but ${dependency}_REQUIRED is ON|AUTO
        # disable module compilation
        if(${dependency}_REQUIRED MATCHES "ON|AUTO")
            string(CONCAT s
                "${message_intro}"
                "Since ${module}_REQUIRED has precedence over ${dependency}_REQUIRED,\ "
                "${module} compilation is disabled.\n"
                "If you do wish to compile ${module}, please set ${module}_REQUIRED to ON or AUTO.
                "
                )
            message(WARNING "${s}")
        endif()
        set(${module}_REQUIRED OFF CACHE STRING "" FORCE)
        set(${module}_FOUND FALSE CACHE INTERNAL "")

    elseif(${module}_REQUIRED STREQUAL "ON")

        # if ${module}_REQUIRED is ON and ${dependency}_REQUIRED is ON:
        #   there is nothing to do

        # if ${module}_REQUIRED is ON but ${dependency}_REQUIRED is AUTO|OFF
        # force ${dependency} to be searched
        if(${dependency}_REQUIRED MATCHES "AUTO|OFF")
            string(CONCAT s
                "${message_intro}"
                "${dependency} is required to compile ${module} and\ "
                "since ${module}_REQUIRED has precedence over ${dependency}_REQUIRED,\ "
                "${dependency}_REQUIRED was changed to ON\ "
                "to enable ${dependency} detection and therefore ${module} compilation.
                "
                )
            message(WARNING "${s}")
        endif()
        set(${dependency}_REQUIRED ON CACHE STRING "" FORCE)

    elseif(${module}_REQUIRED STREQUAL "AUTO")

        # if ${module}_REQUIRED is AUTO and ${dependency}_REQUIRED is ON:
        #   there is nothing to do

        # if ${module}_REQUIRED is AUTO and ${dependency}_REQUIRED is AUTO
        # force ${dependency} to be searched
        if(${dependency}_REQUIRED STREQUAL "AUTO")
            string(CONCAT s
                "${message_intro}"
                "Since ${dependency} is required to compile ${module},\ "
                "${dependency}_REQUIRED was changed to ON\n"
                "to enable ${dependency} detection and therefore ${module} compilation.
                "
                )
            message(STATUS "${s}")
            set(${dependency}_REQUIRED ON CACHE STRING "" FORCE)

        # if ${module}_REQUIRED is AUTO but ${dependency}_REQUIRED is OFF
        # disable module compilation
        elseif(${dependency}_REQUIRED STREQUAL "OFF")
            string(CONCAT s
                "${message_intro}"
                "Since ${dependency} is required to compile ${module},\ "
                "${module}_REQUIRED was changed to OFF to disable ${module} compilation.\n"
                "If you do wish to compile ${module}, please set ${dependency}_REQUIRED to ON.
                "
                )
            message(WARNING "${s}")
            set(${module}_REQUIRED OFF CACHE STRING "" FORCE)
            set(${module}_FOUND FALSE CACHE INTERNAL "")
        endif()

    endif()

    if("${${module}_REQUIRED}" MATCHES "OFF")
        set(${module}_FOUND_MESSAGE
            "Automatically set ${module}_REQUIRED to ${${module}_REQUIRED}"
        CACHE INTERNAL ""
        )
    endif()

endfunction()


function(helyx_set_compiler_flag_for_file file flag)
    if(NOT ${ARGC} EQUAL 2)
        string(CONCAT s
            "helyx_set_compiler_flag_for_file() requires exactly two "
            "arguments, but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    endif()
    set(prev_quiet ${CMAKE_REQUIRED_QUIET})
    set(CMAKE_REQUIRED_QUIET TRUE)
    check_cxx_compiler_flag(${flag} supports${flag})
    if(supports${flag})
        message(STATUS
            "Setting \"${flag}\" for \"${file}\" in target \"${TARGET_NAME}\""
        )
        get_source_file_property(existing_compile_flags ${file}
            COMPILE_FLAGS
        )
        set_source_files_properties(${file}
            PROPERTIES
                COMPILE_FLAGS "${existing_compile_flags} -Wno-unused-private-field"
        )
    endif()
    set(CMAKE_REQUIRED_QUIET ${prev_quiet})
endfunction()



# ============================================================================ #
# ----------------- Functions for configuring HELYX targets ------------------ #
# ============================================================================ #


# ---------------------------- Backend functions ----------------------------- #

# The following targets are intended for internal use - it's not expected that
# users will have to touch these very often

if ("${CMAKE_VERSION}" VERSION_LESS 3.4)
    function(add_library)
        _add_library(${ARGV})
        set(_target_name "${ARGV0}")
        get_target_property(target_type "${_target_name}" TYPE)
        if (NOT "${target_type}" STREQUAL "INTERFACE_LIBRARY")
            set_target_properties(${_target_name}
                PROPERTIES SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
            )
        endif()

        set_directory_properties(PROPERTIES CONTAINS_TARGET_DEFINITION TRUE)
    endfunction()
else()
    function(add_library)
        _add_library(${ARGV})
        set_directory_properties(PROPERTIES CONTAINS_TARGET_DEFINITION TRUE)
    endfunction()
endif()

function(add_executable)
    _add_executable(${ARGV})
    set_directory_properties(PROPERTIES CONTAINS_TARGET_DEFINITION TRUE)
endfunction()

function(add_subdirectory)
    _add_subdirectory(${ARGV})
    set_directory_properties(PROPERTIES ADDS_SUBDIRECTORIES TRUE)
endfunction()


# Directory target names are relative paths from HELYX_PROJECT_DIR to dir, with
# "/" replaced with "-", and -dir appended to the end.
function (get_directory_target_name output_var dir_name)

    # Accept both relative and absolute paths

    # First, check if this is a subdir of the emake starting location
    # Could use OLDPWD env var here, but this seems more in-keeping with the
    # idea that we should avoid env vars where possible
    if (EXISTS "${INITIAL_EMAKE_DIR}/${dir_name}")
        # Get directory as an absolute path (probably unneccesary...)
        get_filename_component(dir_name "${INITIAL_EMAKE_DIR}/${dir_name}" ABSOLUTE)
        # Make relative to project root (may be an external project!)

    endif()

    # Try and resolve this as a relative path from top-level
    # This is probably a bad idea - the intuitive thing is just that users
    # point emake to a path they can actually see, rather than guessing that
    # they might be talking about subdirs of HELYX_PROJECT_DIR
    # if (NOT EXISTS "${CMAKE_SOURCE_DIR}/${dir_name}")
    #     # Get directory as an absolute path (probably unneccesary...)
    #     get_filename_component(ABS_DIR "${CMAKE_SOURCE_DIR}/${dir_name}" ABSOLUTE)
    #     # Make relative to project root (may be an external project!)
    #     file(RELATIVE_PATH dir_name ${CMAKE_SOURCE_DIR} ${ABS_DIR})
    # endif()

    if(IS_ABSOLUTE "${dir_name}")
        file(RELATIVE_PATH dir_name ${CMAKE_SOURCE_DIR} ${dir_name})
    endif()

    string(REPLACE "/" "-" DIR_TARGET_NAME "${dir_name}")
    set(${output_var} "${DIR_TARGET_NAME}-dir" PARENT_SCOPE)
endfunction()


# Assumes that all targets have been added as dependencies of the given dir!
function (get_targets_in_directory output_var dir)
    if("${CMAKE_VERSION}" VERSION_LESS 3.8)
        string(CONCAT s
            "The get_targets_in_directory function is only available in CMake "
            "version 3.8 and above.  The current CMake version is "
            "${CMAKE_VERSION}, so this function will return an empty string."
        )
        message(WARNING "${s}")
        return()
    endif()

    get_directory_target_name(dir_target "${dir}")
    if (NOT TARGET "${dir_target}")
        string(CONCAT s
            "get_targets_in_directory received the following dir, which does "
            "not have a target associated with it:\n"
            "\t${dir}"
            "\t(expecting a target called ${dir_target})"
        )
        message(WARNING "${s}")
        return()
    endif()

    # MANUALLY_ADDED_DEPENDENCIES only available for CMake > 3.8
    get_target_property(temp "${dir_target}" MANUALLY_ADDED_DEPENDENCIES)

    set(${output_var} "${temp}" PARENT_SCOPE)
endfunction()


# Add this target as a dependency of the current directory target, and the
# current directory target as a dependency to all parent directory targets
function (add_to_dir_target TARGET_NAME)
    # Get directory as an absolute path (probably unneccesary...)
    get_filename_component(ABS_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
    # Make relative to project root
    file(RELATIVE_PATH REL_DIR ${CMAKE_SOURCE_DIR} ${ABS_DIR})

    # Add to all higher directories
    string(LENGTH "${REL_DIR}" STR_LEN)
    set(LAST_SLASH_POS ${STR_LEN})
    while (NOT ${LAST_SLASH_POS} LESS 0 AND REL_DIR)
        # Remove after last slash
        string(SUBSTRING ${REL_DIR} 0 ${LAST_SLASH_POS} REL_DIR)
        if (REL_DIR)
            get_directory_target_name(DIR_TARGET_NAME "${REL_DIR}")
            if (NOT TARGET ${DIR_TARGET_NAME})
                add_custom_target(${DIR_TARGET_NAME})

                # Append to list of Helyx targets for emake autocompletion
                # Don't auto-complete directory targets (they're an
                # implementation detail, and there are about 600 of them)
                # get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_DIR_TARGETS)
                # gst(APPEND TARGETS ${DIR_TARGET_NAME})
                # get_property(GLOBAL PROPERTY ALL_HELYX_DIR_TARGETS ${TARGETS})

            endif ()
            add_dependencies(${DIR_TARGET_NAME} ${TARGET_NAME})
            string(FIND ${REL_DIR} "/" LAST_SLASH_POS REVERSE)
        endif ()
    endwhile ()
endfunction ()


macro(add_target_install_logic target_name target_type)
    # Can't add out-of-dir install instructions for CMake < 3.13

    # For Windows, only pack runtime (.exe and .dll, NOT .dll.a)
    if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
        set(RUNTIME_TOGGLE RUNTIME)
    endif()

    if(
        (NOT DEFINED module AND "${module}" STREQUAL "") OR
        ("${module}" MATCHES "^swak4Foam$|^waves2Foam$|^runTimePostProcessing$")
    )
        set(component_name ${HELYX_PACKAGE_NAME})
        set(export_package_name "HELYX-core")
    elseif(NOT "${module}" STREQUAL "")  # All other modules
        set(export_package_name "${module}")
        # In theory, any customer may add any custom module, so don't check names
        set(component_name ${HELYX_PACKAGE_NAME}-${module})
    else()
        string(CONCAT s
            "Failed to add install logic to target \"${target_name}\"\n"
            "The \"module\" keyword was defined (indicating that this target "
            "is not a core target), but it takes the value \"${module}\"."
        )
        message(FATAL_ERROR "${s}")
    endif()

    install(TARGETS ${target_name}
        EXPORT ${export_package_name}
        ${RUNTIME_TOGGLE}
            DESTINATION ${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}/platforms/${HELYX_OPTIONS}/${target_type}
            COMPONENT ${component_name}
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
endmacro()


macro(target_replace_patch_in_exe_name)
    strip_exe_from_target_name(${TARGET_NAME} EXE)
    string(REPLACE "Patch" "Ptch" NEW_NAME ${EXE})
    string(REPLACE "patch" "ptch" NEW_NAME ${NEW_NAME})
    string(COMPARE NOTEQUAL ${EXE} ${NEW_NAME} STRING_CHANGED)
    if (STRING_CHANGED)
        set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${NEW_NAME})
        # Write the batch file at build-time
        add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -Dfilename=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXE}.bat -Dcontent='@${NEW_NAME}.exe %*' -P ${CMAKE_SOURCE_DIR}/etc/cmake/scripts/buildTimeFileWriter.cmake
        )
        # Install the batch file
        if((NOT DEFINED module AND "${module}" STREQUAL "") OR ("${module}" MATCHES "^swak4Foam$|^waves2Foam$|^runTimePostProcessing$|^hisa$"))
            set(component_name ${HELYX_PACKAGE_NAME})
        else()
            set(component_name ${HELYX_PACKAGE_NAME}-${module})
        endif()

        install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXE}.bat
    #         EXPORT export_name  # Will be useful later
            DESTINATION ${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}/platforms/${HELYX_OPTIONS}/bin
            COMPONENT ${component_name}
        )
    endif ()
endmacro()


function(process_windows_executable TARGET_NAME)
    get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
    if ("${TARGET_TYPE}" STREQUAL "EXECUTABLE")
        get_target_property(LINK_LIBRARIES_STRING ${TARGET_NAME} LINK_LIBRARIES)
        set(LIBS_TO_LOAD "")
        foreach(lib ${LINK_LIBRARIES_STRING})
            if(TARGET "${lib}")
                get_target_property(link_target_type "${lib}" TYPE)
                if (link_target_type STREQUAL "SHARED_LIBRARY")
                    list(APPEND LIBS_TO_LOAD "${lib}")
                endif ()
            elseif(EXISTS "${lib}")
                # Probably a library from ThirdParty, e.g.
                # /path/to/libstack_track.dll
                get_filename_component(lib "${lib}" NAME)
                get_filename_component(lib_suffix "${lib}" EXT)
                # Also need to filter for shared libraries here
                if("${lib_suffix}" STREQUAL ${CMAKE_SHARED_LIBRARY_SUFFIX})
                    list(APPEND LIBS_TO_LOAD "${lib}")
                endif()
            else()
                string(CONCAT w
                    "Failed to parse the following for run-time library "
                    "loading (not a target or an existing file):\n\t${lib}"
                )
                message(WARNING "${w}")
            endif ()
        endforeach()
        list(REMOVE_DUPLICATES LIBS_TO_LOAD)
        string(REPLACE ";" "," LIBS_TO_LOAD "${LIBS_TO_LOAD}")
        target_compile_definitions(${TARGET_NAME}
            PRIVATE "LIBS_TO_LOAD=\"${LIBS_TO_LOAD}\""
        )
        unset(LIBS_TO_LOAD)
    endif ()
endfunction()



# -------------------------- User-facing functions --------------------------- #

# It is expected that these functions will be used when a user adds a HELXY
# binary.

function(helyx_link_libraries tgt)
    # Default to private
    set(interface_specifier PRIVATE)
    foreach(linkee ${ARGN})
        if (${linkee} MATCHES "PUBLIC|PRIVATE|INTERFACE")
            set(interface_specifier "${linkee}")
            continue()
        endif()

        target_link_libraries(${tgt}
            ${interface_specifier} ${linkee}
        )

        if (TARGET ${tgt}_obj)
            # This is a silly heuristic to detect file paths vs. targets. Since
            # targets may be created after we have to link to them, we have to
            # decide now, and can't use `if (TARGET ${linkee})`.
            #
            # Things could be neater if the finder modules consistently produced
            # IMPORTED targets instead of a bunch of variables with paths and
            # such, but that proved challenging to make work on windows with the
            # thirdparty stuff.
            if (NOT ${linkee} MATCHES "[/\\]")
                target_include_directories(${tgt}_obj
                    ${interface_specifier} $<TARGET_PROPERTY:${linkee},INTERFACE_INCLUDE_DIRECTORIES>
                )
            endif()
        endif()
    endforeach()
endfunction()

function(helyx_link_libraries_with_generated_headers tgt)
    helyx_link_libraries(${tgt} ${ARGN})

    # Default to private
    set(interface_specifier PRIVATE)
    foreach(linkee ${ARGN})
        if (${linkee} MATCHES "PUBLIC|PRIVATE|INTERFACE")
            set(interface_specifier "${linkee}")
            continue()
        endif()

        if (TARGET ${tgt}_obj)
            target_link_libraries(${tgt}_obj PRIVATE ${linkee})
        endif()
    endforeach()
endfunction()


function(helyx_additional_includes tgt)
    # Default to private
    set(interface_specifier PRIVATE)
    foreach(include ${ARGN})
        if (${include} MATCHES "PUBLIC|PRIVATE|INTERFACE")
            set(interface_specifier "${include}")
            continue()
        endif()

        if("${include}" MATCHES "-includes$")
            target_link_libraries(${tgt}
                "${interface_specifier}" "${include}"
            )

            if (TARGET ${tgt}_obj)
                target_link_libraries(${tgt}_obj
                    "${interface_specifier}" "${include}"
                )
            endif()

            # Keep a note of all un-resolved include targets.  Items from this
            # list are removed when the include target is created.  Remaining
            # items in the list at the end of configuration raise errors.
            if(NOT TARGET "${include}")
                set(UNRESOLVED_LINKS_TO_INCLUDES
                    "${UNRESOLVED_LINKS_TO_INCLUDES};${include}"
                    CACHE INTERNAL
                    ""
                )
                set(test "${${include}-included-from-here}")
                list(APPEND test "${CMAKE_CURRENT_LIST_FILE}")
                set("${include}-included-from-here" "${test}"
                    CACHE INTERNAL ""
                )
            endif()

            continue()
        elseif(TARGET "${include}")
            string(CONCAT s
                "Non-include target linked in add_helyx_library()"
            )
            message(WARNING "${s}")
            target_link_libraries("${tgt}"
                "${interface_specifier}" "${include}"
            )
        endif()

        # Need to get full path for things like IS_DIRECTORY to be well-defined.
        if(NOT IS_ABSOLUTE "${include}")
            get_filename_component(include "${include}" ABSOLUTE)
        endif()

        # It's gotta be a directory or a generator expression, then. Let's
        # hope target_include_directories knows what to do with it.
        target_include_directories("${tgt}"
            "${interface_specifier}" "${include}"
        )

        if (TARGET ${tgt}_obj)
            target_include_directories(${tgt}_obj
                "${interface_specifier}" "${include}"
            )
        endif()
    endforeach()
endfunction()


# TODO:  Some of this logic appears awkwardly placed - perhaps some should be in
# add_helyx_exe and add_helyx_library instead?  Requirement for empty calls to
# this function is very ugly.
# TODO:  Now that lnIncludes are libraries, this logic could be simplified.  It
# could look pretty much identical to the -includes targets.  However, we know
# this works and is well tested, so the best thing to do is instead to just
# deprecate and remove lnIncludes.
#
# This function will become obsolete after the lnIncludes are completely deprecated.
# At that point, user will need to use the helyx_additional_includes() instead.
function(helyx_include_directories TARGET_NAME)

    # By default, things should be PRIVATE
    set(interface_specifier PRIVATE)
    foreach(INCLUDE_DIR ${ARGN})
        # TODO:  Suspect there's a bug here - how would the following (perfectly
        # normal CMake) code be dealt with?
        # helyx_include_directories(${TARGET_NAME}
        #     PRIVATE
        #         thing1
        #         thing2
        #     PUBLIC
        #         thing3
        #         thing4
        # )
        if (${INCLUDE_DIR} MATCHES "PUBLIC|PRIVATE|INTERFACE")
            set(interface_specifier "${INCLUDE_DIR}")
            continue()
        endif()

        # Useful for headers included in a library's header files
        # Add non-lnIncludes as private include directories
        get_filename_component(INCLUDE_DIR ${INCLUDE_DIR} ABSOLUTE)
        file(RELATIVE_PATH RELATIVE_INC_DIR ${HELYX_PROJECT_DIR} ${INCLUDE_DIR})
        target_include_directories(${TARGET_NAME}
            ${interface_specifier}  # One of PUBLIC, INTERFACE, PRIVATE
                $<BUILD_INTERFACE:${HELYX_PROJECT_DIR}/${RELATIVE_INC_DIR}>
                $<INSTALL_INTERFACE:${HELYX_DIR_NAME}/${RELATIVE_INC_DIR}>
        )

        # Reset the interface specifier (N.B. to get here, the continue
        # statement must have been missed)
        set(interface_specifier PRIVATE)
    endforeach()

    string(CONCAT s
    "LnIncludes folders have been deprecated in HELYX. "
    "This message is intended principally for custom codes compilation against HELYX. \n"
    "HELYX-Core and Modules are already independent of lnIncludes but, "
    "in the current version, custom codes can still be compiled against HELYX using the lnIncludes. "
    "With the complete deprecation of lnIncludes in the next releases, "
    "the helyx_include_directory() function will become obsolete "
    "and users will need to update their custom codes accordingly.\n"
    "More details about compiling custom codes against HELYX without lnIncludes "
    "can be found in the 'HELYX-Core Compilation Guide'.\n"
    )
    message(AUTHOR_WARNING "${s}")

endfunction()


function(add_helyx_include_target top_level_target_dir)
    # Need relative paths for installation
    get_filename_component(relative_inc_dir ${top_level_target_dir} ABSOLUTE)
    file(RELATIVE_PATH relative_inc_dir
        ${HELYX_PROJECT_DIR} ${relative_inc_dir}
    )

    add_library(${TARGET_NAME}-includes INTERFACE)
    target_include_directories(${TARGET_NAME}-includes BEFORE
        INTERFACE
            $<BUILD_INTERFACE:${top_level_target_dir}>
            $<INSTALL_INTERFACE:${HELYX_DIR_NAME}/${relative_inc_dir}>
    )
    get_target_property(target_type "${TARGET_NAME}" TYPE)
    if (target_type MATCHES "^INTERFACE_LIBRARY$")
        target_link_libraries(${TARGET_NAME} INTERFACE ${TARGET_NAME}-includes)
    else()
        target_link_libraries(${TARGET_NAME} PUBLIC ${TARGET_NAME}-includes)
    endif ()

    add_target_install_logic("${TARGET_NAME}-includes" "include_targets")

    set(unresolved_includes "${UNRESOLVED_LINKS_TO_INCLUDES}")
    list(REMOVE_DUPLICATES unresolved_includes)
    list(REMOVE_ITEM unresolved_includes "${TARGET_NAME}-includes")
    set(UNRESOLVED_LINKS_TO_INCLUDES
        "${unresolved_includes}"
        CACHE INTERNAL
        ""
    )
endfunction()


# TODO:  Consider moving some of this functionality to an add_library() override
# so targets like OpenFOAM don't need to duplicate quite so much code
function(add_helyx_library TARGET_NAME TARGET_TYPE)
    add_library(${TARGET_NAME}_obj OBJECT ${ARGN})

    add_library(${TARGET_NAME} ${TARGET_TYPE})
    target_link_libraries(${TARGET_NAME} PRIVATE ${TARGET_NAME}_obj)

    target_include_directories(${TARGET_NAME}_obj
        PRIVATE $<TARGET_PROPERTY:${TARGET_NAME},INCLUDE_DIRECTORIES>
    )

    # Although strictly speaking this is handled by the includes target, in
    # order to specify that the PWD should take precendence over everything
    # else, we have to specify this here with the BEFORE keyword.  This appears
    # not to be necessary for libraries, but better safe than sorry...
    # We still require the includes targets for header-only code.
    target_include_directories(${TARGET_NAME}_obj BEFORE
        PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    )

    # Append to list of Helyx targets for emake autocompletion
    get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
    list(APPEND TARGETS ${TARGET_NAME})
    set_property(GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS ${TARGETS})

    add_helyx_include_target(${CMAKE_CURRENT_LIST_DIR})

    add_to_dir_target(${TARGET_NAME})

    if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
        if (DEFINED module AND NOT "${module}" STREQUAL "")
            set(OUTPUT_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${module}")
        else()
            set(OUTPUT_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
        endif()
        # CMake treats Windows .dll files as runtime targets
        set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${OUTPUT_DIR}")
    endif()

    # Only pack shared targets
    if(${TARGET_TYPE} STREQUAL "SHARED")
        add_target_install_logic("${TARGET_NAME}" lib)
    endif()

    if (${HELYX_SYSTEM_NAME} STREQUAL "POSIX")
        target_link_options(${TARGET_NAME}
            PRIVATE -ldl -lm
        )
    endif()

    get_target_property(install_runpath ${TARGET_NAME} INSTALL_RPATH)
    string(REPLACE "$ORIGIN/../lib"
        "$ORIGIN" install_runpath
        "${install_runpath}"
    )
    set_target_properties(${TARGET_NAME}
        PROPERTIES INSTALL_RPATH "${install_runpath}"
    )

endfunction()


function(add_helyx_executable TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})

    # Although strictly speaking this is handled by the includes target, in
    # order to specify that the PWD should take precendence over everything
    # else, we have to specify this here with the BEFORE keyword.  This is
    # required for exe's, because exe's often use horrible things like
    # "createFields.H", and if the order of includes is incorrect then things
    # can go horribly wrong.
    # We still require the includes targets for header-only code.
    target_include_directories(${TARGET_NAME} BEFORE
        PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    )

    add_helyx_include_target(${CMAKE_CURRENT_LIST_DIR})

    # Append to list of Helyx targets for Windows runtime library loading and
    # emake autocompletion
    get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_EXE_TARGETS)
    list(APPEND TARGETS ${TARGET_NAME})
    set_property(GLOBAL PROPERTY ALL_HELYX_EXE_TARGETS ${TARGETS})

    add_target_install_logic("${TARGET_NAME}" bin)

    # This is now handled by the include target.
    # target_include_directories(${TARGET_NAME}
    #     PRIVATE .
    # )

    add_to_dir_target(${TARGET_NAME})

    # Override output directory for modules
    if (${HELYX_SYSTEM_NAME} STREQUAL "MSwindows" AND DEFINED module AND NOT "${module}" STREQUAL "")
        set(OUTPUT_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${module}")
    endif()

    strip_exe_from_target_name(${TARGET_NAME} EXE_NAME)
    # Set exe name
    set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${EXE_NAME})

    if ("${HELYX_SYSTEM_NAME}" STREQUAL "POSIX")
        target_link_options(${TARGET_NAME}
            PRIVATE -ldl -lm
        )
    elseif("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
        target_replace_patch_in_exe_name()
    endif ()
endfunction()


function(pad_string_with_char string_to_pad character max_length)

    # If centred text, stop one short rather than overunning
    if("${ARGN}" STREQUAL "centre")
        math(EXPR max_length ${max_length}-1)
    endif()

    set(string_to_print "${string_to_pad}")
    string(LENGTH ${string_to_print} string_to_print_length)
    while(${string_to_print_length} LESS ${max_length})
        if("${ARGN}" STREQUAL "centre")
            set(string_to_print "${character}${string_to_print}")
        endif()
        set(string_to_print "${string_to_print}${character}")
        string(LENGTH ${string_to_print} string_to_print_length)
    endwhile()

    if("${ARGN}" STREQUAL "centre")
        math(EXPR max_length ${max_length}+1)
    endif()

    # This is required for centred text, and harmless otherwise
    if(${string_to_print_length} LESS ${max_length})
        set(string_to_print "${string_to_print}${character}")
    endif()

    set(string_to_print "${string_to_print}" PARENT_SCOPE)
endfunction()


function(helyx_flex_target target_name)
    flex_target(${target_name} ${ARGN})
    # Reduce warning level for flex generated file
    set(OUTPUTS FLEX_${target_name}_OUTPUTS)
    # In newer versions of cmake, COMPILE_OPTIONS can be used to set a list
    # (instead of space-separated flags currently used)
    set_source_files_properties(${${OUTPUTS}} PROPERTIES COMPILE_FLAGS "${LESSWARN}")
    set(${OUTPUTS} ${${OUTPUTS}} PARENT_SCOPE)
endfunction()






function(get_all_source_files output_var top_level_dir)
    # message("top_level_dir = ${top_level_dir}")
    set(include_dirs "${top_level_dir}")

    # Match files (this is probably daft, but backwards-compatability with
    # lnIncludes seems like a good idea at this stage).
    file(GLOB_RECURSE files_in_subdirs
        LIST_DIRECTORIES FALSE
        # RELATIVE ${HELYX_PROJECT_DIR}
        ${top_level_dir}/*.[CHhL]
        ${top_level_dir}/*.[ch]xx
        ${top_level_dir}/*.[ch]pp
        ${top_level_dir}/*.type
    )

    set(output "")

    foreach(subdir_file ${files_in_subdirs})
        string(REGEX MATCH lnInclude\/
            test ${subdir_file}
        )
        # message("subdir_file = ${subdir_file}")
        if(NOT "" STREQUAL "${test}")
            list(REMOVE_ITEM files_in_subdirs ${subdir_file} )
        endif()
        if("" STREQUAL "${test}")
            get_filename_component(subdir_file ${subdir_file} REALPATH)
            list(APPEND output ${subdir_file})
        endif()
    endforeach()
    list(REMOVE_DUPLICATES output)

    # TODO:  Think about duplicates...
    # list(APPEND list_of_files "${files_in_subdirs}")

    # set(${output_var} "${include_dirs}" PARENT_SCOPE)
    set(${output_var} "${output}" PARENT_SCOPE)
endfunction()


function(get_includes_in_file output_var input_file)
    file(READ ${input_file} input_string)
    string(REGEX MATCHALL
        # Assume system includes use <>, OF includes use ""
        "#[ \t\r\n]*include[ \t\r\n]+\"([^\"\n]*)\""
        all_matches
        "${input_string}"
    )
    # CMake seems not to be setting the group matches correctly, so just treat
    # this as a "#include "-separated list

    string(REGEX REPLACE
        "#[ \t\r\n]*include[ \t\r\n]+\"([^\"\n]*)\""
        "\\1"
        all_matches
        "${all_matches}"
    )
    list(REMOVE_ITEM all_matches
        "scotch.h"
        "mpi.h"
    )
    set(${output_var} "${all_matches}" PARENT_SCOPE)
endfunction()


function(resolve_include output_var original_source_file file_to_be_included include_dirs)
    if(NOT ${ARGC} EQUAL 4)
        string(CONCAT s
            "resolve_include() requires exactly four arguments, but ${ARGC} "
            "were provided"
        )
        message(FATAL_ERROR "${s}")
    endif()

    # message("${original_source_file}:  ${file_to_be_included}")

    # Sometime the include is in the same dir as the source
    get_filename_component(src_file_dir "${original_source_file}" DIRECTORY)
    if(EXISTS "${src_file_dir}/${file_to_be_included}")
        set(${output_var} "${file_to_be_included}" PARENT_SCOPE)
        # message("    FOUND ${file_to_be_included} in same dir")
        # message("    Path  ${src_file_dir}")
        return()
    endif()

    # This is ugly, but we have to deal with generator expressions...  Could
    # also consider making targets to do this?  But would be harder to debug...
    string(REGEX REPLACE "\\\$<BUILD_INTERFACE:([^>]*)>" "\\1"
        include_dirs "${include_dirs}")
    string(REGEX REPLACE "\\\$<INSTALL_INTERFACE:([^>]*)>" ""
        include_dirs "${include_dirs}")

    foreach(search_dir ${include_dirs})
        # For now, ignore lnInclude directories - they exist in includes_dirs,
        # so we should either not include them in there or we should actually
        # use them
        string(REGEX MATCH "lnInclude/?$"
            match "${search_dir}")
        if(NOT "" STREQUAL "${match}")
            # message("Ignoring lnInclude dir ${search_dir}")
            continue()
        endif()


        if(NOT IS_ABSOLUTE "${search_dir}")
            # Make absolute for consistency
            get_filename_component(search_dir "${search_dir}"
                ABSOLUTE BASE_DIR "${src_file_dir}"
            )
        endif()

        if(NOT EXISTS "${search_dir}")
            # Could be an lnInclude that hasn't been generated
            continue()
        else()
            # Don't mess with paths outside HELYX_PROJECT_DIR
            get_filename_component(real_search_dir "${search_dir}" REALPATH)
            get_filename_component(real_proj_dir "${HELYX_PROJECT_DIR}" REALPATH)
            if(NOT "${real_search_dir}" MATCHES "${HELYX_PROJECT_DIR}")
                continue()
            endif()
        endif()

        # find_file isn't terribly helpful, because it doesn't recurse
        # unset(full_path_to_file)
        # find_file (full_path_to_file
        #     "${file_to_be_included}"
        #     HINTS "${search_dir}"
        #     NO_DEFAULT_PATH
        #     NO_CACHE
        # )

        if(EXISTS "${search_dir}/${file_to_be_included}")
            # This is already a valid include!  Do nothing
            set(${output_var} "${file_to_be_included}" PARENT_SCOPE)
            return()
        else()
            list(APPEND USED_SEARCH_DIRS "${search_dir}")
            set(USED_SEARCH_DIRS "${USED_SEARCH_DIRS}" PARENT_SCOPE)
            # file(GLOB_RECURSE full_path_to_file
            file(GLOB_RECURSE temp
                RELATIVE "${search_dir}"
                FOLLOW_SYMLINKS
                LIST_DIRECTORIES FALSE
                "${search_dir}/*/${file_to_be_included}"
            )
            set(full_path_to_file "")
            foreach(item ${temp})
                if(NOT "${item}" MATCHES lnInclude)
                    list(APPEND full_path_to_file "${item}")
                endif()
            endforeach()

            list(LENGTH full_path_to_file len)
            if(1 LESS "${len}")
                string(CONCAT s
                    "Glob picked up more than one file!\n"
                    "    File name:  \"${file_to_be_included}\"\n"
                    "    Glob pattern:  \"${search_dir}/*/${file_to_be_included}\""
                    "    Result:  \"${full_path_to_file}\""
                )
                message(SEND_ERROR "${s}")
            endif()

            if(NOT "" STREQUAL "${full_path_to_file}"
                AND EXISTS "${search_dir}/${full_path_to_file}"
            )
                set(${output_var} "${full_path_to_file}" PARENT_SCOPE)
                return()
            endif()
        endif()
    endforeach()

    # message("Failed to find include!\n\n")
    set(${output_var} "" PARENT_SCOPE)
endfunction()


function(replace_include source_file original_include full_path)
    # Shouldn't be possible to get in here anymore
    if(NOT full_path)
        # Include may be in the same dir as the source file
        get_filename_component(src_file_dir "${source_file}" DIRECTORY)
        if(NOT EXISTS "${src_file_dir}/${original_include}")
            string(CONCAT s
                "Got the following invalid refactoring information, not doing "
                "anything:\n"
                "\tsource_file:\t${source_file}\n"
                "\toriginal_include:\t${original_include}\n"
                "\tfull_path:\t${full_path}"
            )
               message(FATAL_ERROR "${s}")
            return()
        endif()
    endif()

    file(READ ${source_file} input_string)

    string(REGEX REPLACE
        "#[ \t\r\n]*include[ \t\r\n]+\"${original_include}\""
        "#include \"${full_path}\""
        input_string
        "${input_string}"
    )

    file(WRITE ${source_file} "${input_string}")
endfunction()


function(refactor_includes include_dirs)
    get_all_source_files(source_files "${CMAKE_CURRENT_LIST_DIR}")

    foreach(source_file ${source_files})
        get_includes_in_file(includes ${source_file})

        foreach(include ${includes})
            resolve_include(full_path_to_include "${source_file}" "${include}" "${include_dirs}")
            if(NOT "${full_path_to_include}" STREQUAL "")
                replace_include("${source_file}" "${include}" "${full_path_to_include}")
            else()
                message(WARNING "full_path_to_include was empty
                source_file: \"${source_file}\"
                include: \"${include}\"
                full_path_to_include: \"${full_path_to_include}\"
                "
                )
            endif()
        endforeach()

    endforeach()
endfunction()


# Could definitely do something nicer here (e.g. memoization), but that seems
# like a premature optimisation.
#
# Bug: There is CMake bug for version lower than 3.19.5 when the TARGET_NAME is an INTERFACE_LIBRARY
#   INTERFACE_LIBRARY targets may only have whitelisted properties.  The
#   property "LINK_LIBRARIES" is not allowed.
#
# Since this function is used only for the Cpp and the CMakeLists inludes refactor and lnIncludes
# will be completere deprecated in the next releases, for now the solution is simply to warn the user to use
# a CMake version >= 3.19.5. If need another fix in the future we can use
# get_target_property(target_type "${TARGET_NAME}" TYPE) to skip the INTERFACE_LIBRARY targets
function(target_get_all_target_dependencies output_var TARGET_NAME)
    set(dependencies "")
    # message("Getting direct dependencies of ${TARGET_NAME}...")

    get_target_property(linked_libs "${TARGET_NAME}" LINK_LIBRARIES)
    foreach(linked_lib ${linked_libs})
        if(TARGET "${linked_lib}")
            # message("    ${linked_lib}")
            list(APPEND dependencies "${linked_lib}")
            target_get_all_target_dependencies(next_level_dependencies "${linked_lib}")
            if(NOT "" STREQUAL "${next_level_dependencies}")
                list(APPEND dependencies "${next_level_dependencies}")
            endif()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES dependencies)
    set(${output_var} "${dependencies}" PARENT_SCOPE)
endfunction()


function(check_mpicc mpicc_variable_name)

    find_program (MPICC
        NAMES mpicc
        HINTS ${MPI_ARCH_PATH} ${MPI_LIBRARY_PATHS}
        PATH_SUFFIXES bin;../bin
        DOC "Path on which mpiexec is found"
        )

    if (NOT EXISTS "${MPICC}")
        message(SEND_ERROR
            "Failed to find mpicc, your PATH may be set incorrectly. "
            "Please address this before running HELYX."
        )
        return()
    elseif(NOT "${MPICC}" MATCHES "${MPI_ARCH_PATH}")
        message(SEND_ERROR
            "\nThe variable \"${MPICC}\" specifies the following path, "
            "which is not in the \"MPI_ARCH_PATH\": \"${MPICC}\"\n"
            "Please address this before running HELYX."
        )
        return()
    endif()


    # Here follows a very simple compilation check with the specified mpicc.
    # I feel like there is so much that could go wrong with mpicc that these
    # are better off as warnings than errors.

    # Create the mpicc test
    set(output_string "int main () { return 21 + 21\; }")
    if(EXISTS${PROJECT_BINARY_DIR}/testMPICC)
        file(REMOVE_RECURSE ${PROJECT_BINARY_DIR}/testMPICC/)
    endif()
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/testMPICC)
    file(WRITE ${PROJECT_BINARY_DIR}/testMPICC/testMpicc.c ${output_string})

    # Puts .o files in the current binary dir
    execute_process(
        COMMAND ${MPICC} -c testMpicc.c
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/testMPICC
        TIMEOUT 10
        RESULT_VARIABLE test
    )
    if(NOT "0" STREQUAL "${test}")
        message(WARNING
            "Failed to compile a simple program using the following mpicc (specified in the \"${mpicc_variable_name}\" variable):
            \"${MPICC}\"")
        return()
    endif()

    # If we get here, the compilation must have succeeded, so try linking
    execute_process(
        COMMAND ${MPICC} -o testMpicc testMpicc.o
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/testMPICC
        TIMEOUT 10
        RESULT_VARIABLE test
    )
    if(NOT "0" STREQUAL "${test}")
        message(WARNING
            "Failed to link a simple program using the following mpicc (specified in the \"${mpicc_variable_name}\" variable):
            \"${MPICC}\"")
        return()
    endif()

    # If we get here, then we should have an executable...
    execute_process(
        COMMAND ${PROJECT_BINARY_DIR}/testMPICC/testMpicc
        TIMEOUT 10
        RESULT_VARIABLE test
    )
    if(NOT "42" STREQUAL "${test}")
        message(WARNING
            "Got an unexpected result from a simple test program that was compiled with the following mpicc (specified in the \"${mpicc_variable_name}\" variable):
            \"${MPICC}\"")
        return()
    #else()
    #    string(CONCAT s
    #        "The mpicc at the following location compiled a simple test "
    #        "program successfully:\n\t\"${MPICC}\""
    #    )
    #    message(STATUS "${s}")
    endif()

endfunction()

# Add a linker flag, if the linker supports it.
function (optional_linker_flag FLAG)
    if (CMAKE_VERSION VERSION_LESS "3.18")
        # Can't be bothered to figure this out for people who want to use a >5 year old version of cmake.
        message("Upgrade cmake to at least version 3.18 for slightly more optimisation")
        return()
    endif()

    include(CheckLinkerFlag)

    string(MAKE_C_IDENTIFIER ${FLAG} CACHE_VAR)
    check_linker_flag(CXX LINKER:${FLAG} ${CACHE_VAR})
    if (${${CACHE_VAR}})
        add_link_options(LINKER:${FLAG})
    endif()
endfunction()

# Add a compiler flag, if the compiler supports it.
function (optional_compiler_flag _F)
    if (CMAKE_VERSION VERSION_LESS "3.19")
        # Can't be bothered to figure this out for people who want to use a >5 year old version of cmake.
        message("Upgrade cmake to at least version 3.19 for slightly more optimisation")
        return()
    endif()

    include(CheckCompilerFlag)

    string(MAKE_C_IDENTIFIER ${_F} CACHE_VAR)
    check_cxx_compiler_flag(${_F} ${CACHE_VAR})

    if (${${CACHE_VAR}})
        add_compile_options(${_F})
    endif ()
endfunction()

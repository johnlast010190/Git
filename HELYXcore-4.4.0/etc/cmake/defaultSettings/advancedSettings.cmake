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
    (c) 2019-2020 Engys Ltd.

Description
    This file sets defaults for advanced and non-critical core variables, and
    checks that certain critical variables have sensible values.

[----------------------------------------------------------------------------]]


# ============================================================================ #
# --------------------------  Configuration types  --------------------------- #
# ============================================================================ #

set(ENABLE_FPE
    "On"
    CACHE STRING
    "Enable FPE exception"
    FORCE
)

set(CMAKE_CONFIGURATION_TYPES Opt;Prof;Debug;Asan
    CACHE STRING
    "Build types, normally Opt(imised), Prof(ile), Debug"
    FORCE
    )


# Runtime debugging code. Enable by default only for debug builds.
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(_DEFAULT_HELYX_RUNTIME_DEBUG_FLAGS On)
else()
    set(_DEFAULT_HELYX_RUNTIME_DEBUG_FLAGS Off)
endif()

if ("${ENABLE_FPE}" STREQUAL "On")
    set(_DEFAULT_HELYX_RUNTIME_DEBUG_FLAGS On)
endif()

option(HELYX_RUNTIME_DEBUG_FLAGS "Enable runtime debuging code." ${_DEFAULT_HELYX_RUNTIME_DEBUG_FLAGS})

# General settings, applied to all configurations
set(CMAKE_CXX_FLAGS
    "-m64 -ftemplate-depth=1024 -Wall -Wextra"
    CACHE STRING
    "Flags used by the C++ compiler for all builds."
    FORCE
    )
set(CMAKE_CXX_OPTIONAL_FLAGS
    "-Wnon-virtual-dtor -Wno-unused-parameter -Wno-overloaded-virtual -Wno-undefined-var-template -Wno-old-style-cast"
    CACHE STRING
    "Flags used by the C++ compiler, if supported, for all builds"
    FORCE
    )
set(CMAKE_C_FLAGS ""
    CACHE STRING
    "Flags used by the C compiler for all builds."
    FORCE
    )
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed"
    CACHE STRING
    "Flags used for linking binaries for all builds."
    FORCE
    )
set(CMAKE_SHARED_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed"
    CACHE STRING
    "Flags used by the shared libraries linker for all builds."
    FORCE
)

mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS
    CMAKE_C_FLAGS
    CMAKE_EXE_LINKER_FLAGS
    CMAKE_SHARED_LINKER_FLAGS
    )


# Opt-specific settings
set(CMAKE_CXX_FLAGS_OPT "-O3 -DNDEBUG"
    CACHE STRING
    "Flags used by the C++ compiler during \"Opt\" builds."
    FORCE
    )
set(CMAKE_C_FLAGS_OPT "-DNDEBUG"
    CACHE STRING
    "Flags used by the C compiler during \"Opt\" builds."
    FORCE
    )
set(CMAKE_EXE_LINKER_FLAGS_OPT ""
    CACHE STRING
    "Flags used for linking binaries during \"Opt\" builds."
    FORCE
    )
set(CMAKE_SHARED_LINKER_FLAGS_OPT ""
    CACHE STRING
    "Flags used by the shared libraries linker during \"Opt\" builds."
    FORCE
    )
mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS_OPT
    CMAKE_C_FLAGS_OPT
    CMAKE_EXE_LINKER_FLAGS_OPT
    CMAKE_SHARED_LINKER_FLAGS_OPT
    )


# Prof-specific settings
set(CMAKE_CXX_FLAGS_PROF "-pg -O2 -DNDEBUG"
    CACHE STRING
    "Flags used by the C++ compiler during \"Prof\" builds."
    FORCE
    )
set(CMAKE_C_FLAGS_PROF "-DNDEBUG"
    CACHE STRING
    "Flags used by the C compiler during \"Prof\" builds."
    FORCE
    )
set(CMAKE_EXE_LINKER_FLAGS_PROF ""
    CACHE STRING
    "Flags used for linking binaries during \"Prof\" builds."
    FORCE
    )
set(CMAKE_SHARED_LINKER_FLAGS_PROF ""
    CACHE STRING
    "Flags used by the shared libraries linker during \"Prof\" builds."
    FORCE
    )
mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS_PROF
    CMAKE_C_FLAGS_PROF
    CMAKE_EXE_LINKER_FLAGS_PROF
    CMAKE_SHARED_LINKER_FLAGS_PROF
    )


# Debug-specific settings
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0 -fdefault-inline -DFULLDEBUG"
    CACHE STRING
    "Flags used by the C++ compiler during \"Debug\" builds."
    FORCE
    )
set(CMAKE_C_FLAGS_DEBUG ""
    CACHE STRING
    "Flags used by the C compiler during \"Debug\" builds."
    FORCE
    )
set(CMAKE_EXE_LINKER_FLAGS_DEBUG ""
    CACHE STRING
    "Flags used for linking binaries during \"Debug\" builds."
    FORCE
    )
set(CMAKE_SHARED_LINKER_FLAGS_DEBUG ""
    CACHE STRING
    "Flags used by the shared libraries linker during \"Debug\" builds."
    FORCE
    )
mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_C_FLAGS_DEBUG
    CMAKE_EXE_LINKER_FLAGS_DEBUG
    CMAKE_SHARED_LINKER_FLAGS_DEBUG
    )

# Asan-specific settings
set(CMAKE_CXX_FLAGS_ASAN
    "-g -gcolumn-info -O3 -fno-omit-frame-pointer -fsanitize=undefined -fsanitize-address-use-after-scope -fsanitize=address"
    CACHE STRING
    "Flags used by the C++ compiler during \"Asan\" builds."
    FORCE
    )
set(CMAKE_C_FLAGS_ASAN ""
    CACHE STRING
    "Flags used by the C compiler during \"Asan\" builds."
    FORCE
    )
set(CMAKE_EXE_LINKER_FLAGS_ASAN "-fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address"
    CACHE STRING
    "Flags used for linking binaries during \"Asan\" builds."
    FORCE
    )
set(CMAKE_SHARED_LINKER_FLAGS_ASAN "-fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address"
    CACHE STRING
    "Flags used by the shared libraries linker during \"Asan\" builds."
    FORCE
    )
mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS_ASAN
    CMAKE_C_FLAGS_ASAN
    CMAKE_EXE_LINKER_FLAGS_ASAN
    CMAKE_SHARED_LINKER_FLAGS_ASAN
)


# ============================================================================ #
# -----------------------------  CMAKE_COMPILER ------------------------------ #
# ============================================================================ #

# Make sure that CMAKE_CXX_COMPILER and CMAKE_C_COMPILER get the compiler full
# path and not only the name when stored in the CMakeCache.
#
# There is a chance that userSettings.cmake is reloaded
# during auto-configuration, and then the full paths detected by CMake may be
# overriden if only the name is given to CMAKE_*_COMPILER, which will change the
# build rules and, consequently, trigger the whole build!
#
# Need to do it here and not from anything in "projectHeader.cmake"
# because this script is usually called before
# "project(HELYXcore)" (which triggers the CMake system detection)

# We could introduce HELYX_*_COMPILER variables as a  friendly interface to set
# the CC and CXX environment variables, which will in fact be used by CMake to
# detect and populate CMAKE_*_COMPILER variables, as below.
# But this would add extra complexity as we will need to rename variables here
# and in the generateUserSettingsFile.cmake
#if (NOT "" STREQUAL "${HELYX_C_COMPILER}" AND
#    NOT "" STREQUAL "${HELYX_CXX_COMPILER}")
#    # To set the environment variables CC and CXX for CMake system-detection
#    # to not mess up with CMAKE_*_COMPILER variables
#    set(ENV{CXX} "${CMAKE_CXX_COMPILER}")
#    set(ENV{CC} "${CMAKE_C_COMPILER}")

if("" STREQUAL "${CMAKE_C_COMPILER}" OR "" STREQUAL "${CMAKE_CXX_COMPILER}")
    # Probably we will never reach this point
    message(FATAL_ERROR
        "CMAKE_<lang>_COMPILER not found or not defined in your userSettings. "
        "Please set CMAKE_C_COMPILER and CMAKE_CXX_COMPILER and run 'emake -r'\n")

elseif(NOT EXISTS "${CMAKE_C_COMPILER}" OR NOT EXISTS "${CMAKE_CXX_COMPILER}")

    execute_process(COMMAND which "${CMAKE_C_COMPILER}"
        OUTPUT_VARIABLE c_compiler_path
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    execute_process(COMMAND which "${CMAKE_CXX_COMPILER}"
        OUTPUT_VARIABLE cxx_compiler_path
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    # Overwrite the CMAKE_*_COMPILER cache variables with full path
    if(c_compiler_path AND cxx_compiler_path)
        set(CMAKE_C_COMPILER
            "${c_compiler_path}"
            CACHE PATH
            "The C compiler to use (can be a path to an executable)"
            FORCE
            )
        set(CMAKE_CXX_COMPILER
            "${cxx_compiler_path}"
            CACHE PATH
            "The C compiler to use (can be a path to an executable)"
            FORCE
            )
    else()
        message(WARNING "CMake could not resolve the full path of the compiler. "
            "This may lead to conflicts during CMake re-configuration. "
            "Please review the compiler path set on the variables "
            "CMAKE_C_COMPILER and CMAKE_CXX_COMPILER. "
            "If the path is correct, please disregard this message.")
    endif()
endif()


# ============================================================================ #
# ---------------------------  HELYX_COMPILER_NAME --------------------------- #
# ============================================================================ #

# First, set the CXX and CC environment variables, as set out in the CMake FAQ
# set(ENV{CXX} "${CMAKE_CXX_COMPILER}")
# set(ENV{CC} "${CMAKE_C_COMPILER}")

# Can't use CMAKE_C_COMPILER_VERSION, because it's not set when this file is
# called as a script from emake
if ("" STREQUAL "${HELYX_COMPILER_NAME}")
    execute_process(COMMAND
        ${CMAKE_C_COMPILER} -dumpfullversion -dumpversion
        RESULT_VARIABLE RETURN_CODE
        OUTPUT_VARIABLE VERSION_STRING
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT "${RETURN_CODE}" STREQUAL "0")
        message(FATAL_ERROR "Error when setting HELYX_COMPILER_NAME!
        Attempted to get the version from the following C compiler:
            \"${CMAKE_C_COMPILER}\"
        Using the following command:
            \"${CMAKE_C_COMPILER} -dumpfullversion -dumpversion\"
        Received the following error:
            \"${RETURN_CODE}\"
        HELYX_COMPILER_NAME is a purely cosmetic variable (it's only used to set paths),
        but it is necessary.  Please check your C compiler path (CMAKE_C_COMPILER), or
        set HELYX_COMPILER_NAME manually.")
    elseif(VERSION_STRING)
        # Set HELYX_COMPILER_NAME (ONLY useful for setting options correctly)
        get_filename_component(C_COMPILER_NAME "${CMAKE_C_COMPILER}" NAME)
        string(SUBSTRING ${C_COMPILER_NAME} 0 1 FIRST_LETTER)
        string(TOUPPER ${FIRST_LETTER} FIRST_LETTER)
        string(REGEX REPLACE
            "^.(.*)"   "${FIRST_LETTER}\\1"
            CAPITALISED_CMAKE_C_COMPILER
            "${C_COMPILER_NAME}"
            )
        string(REPLACE "." ";" VERSION_LIST "${VERSION_STRING}")
        list(GET VERSION_LIST 0 COMILER_VERSION_MAJOR)
        list(GET VERSION_LIST 1 COMILER_VERSION_MINOR)
        set(HELYX_COMPILER_NAME
            ${CAPITALISED_CMAKE_C_COMPILER}${COMILER_VERSION_MAJOR}${COMILER_VERSION_MINOR}
            CACHE STRING
            "User-friendly name for CMAKE_C_COMPILER (ONLY used for setting directory names)"
        )
    else()
        message(FATAL_ERROR "Unknown error when setting HELYX_COMPILER_NAME!
        Attempted to get the version from the following C compiler:
            \"${CMAKE_C_COMPILER}\"
        Using the following command:
            \"${CMAKE_C_COMPILER} -dumpfullversion -dumpversion\"
        Received the following from stdout:
            \"${VERSION_STRING}\"
        Received the following error:
            \"${RETURN_CODE}\"
        HELYX_COMPILER_NAME is a purely cosmetic variable (it's only used to set paths),
        but it is necessary.  Please set it manually.")
    endif()
endif()


# ============================================================================ #
# ----------------------------------  Paths ---------------------------------- #
# ============================================================================ #

# ------------------------------  HELYX_OPTIONS ------------------------------ #

# This is equivalent to WM_OSTYPE, but renamed to reflect the CMake variable
# "CMAKE_SYSTEM_NAME"
# Consider removing this, and just using CMAKE_SYSTEM_NAME instead
# This should be reset in the toolchain files as necessary
# set(HELYX_SYSTEM_NAME POSIX CACHE STRING "Target system (POSIX or MSwindows)")
set(HELYX_SYSTEM_NAME POSIX
    CACHE STRING
    "Taret system (POSIX or MSwindows)"
    FORCE
    )
mark_as_advanced(FORCE HELYX_TARGET_PLATFORM)

# Set things necessary to set HELYX_OPTIONS
# This is equivalent to WM_ARCH, can be set to linuxmingw_w64 when
# cross-compiling
set(HELYX_BUILD_PLATFORM linux CACHE STRING
    "The platform used for the build (normally linux)"
    FORCE
    )
mark_as_advanced(FORCE HELYX_BUILD_PLATFORM)

# This is confusingly named, but there's no CMake alternative.  We only support
# 64-bit architectures.
# I prefer HELYX_SYSTEM_ADDRESS_SIZE (strictly speaking I think it's closer to
# the OS word size, but that's confusing in an OpenFOAM context...)
# Build-system \"bitness\" (N.B. HELYX only supports 64-bit)
set(HELYX_COMPILER_LIB_ARCH 64
    CACHE INTERNAL
    "The bitness of the build system (HELYX only supports 64-bit systems)"
    FORCE
    )

string(CONCAT TEMP
    "${HELYX_BUILD_PLATFORM}"
    "${HELYX_COMPILER_LIB_ARCH}"
    "${HELYX_COMPILER_NAME}"
    "${HELYX_PRECISION_OPTION}"
    "Int${HELYX_LABEL_SIZE}"
    "${CMAKE_BUILD_TYPE}"
)
set(HELYX_OPTIONS ${TEMP}
    CACHE STRING
    "A string used to set paths in platforms directories (not user-settable)"
    FORCE
    )
mark_as_advanced(HELYX_OPTIONS)

string(CONCAT TEMP
    "${HELYX_BUILD_PLATFORM}"
    "${HELYX_COMPILER_LIB_ARCH}"
    "${HELYX_COMPILER_NAME}"
)
set(HELYX_OPTIONS_COMPILER ${TEMP}
    CACHE STRING
    "A string used to set paths in platforms directories (not user-settable):
    \${HELYX_BUILD_PLATFORM}\${HELYX_COMPILER_LIB_ARCH}\${HELYX_COMPILER_NAME}"
    FORCE
    )
mark_as_advanced(HELYX_OPTIONS_COMPILER)

string(CONCAT TEMP
    "${HELYX_PRECISION_OPTION}"
    "Int${HELYX_LABEL_SIZE}"
)
set(HELYX_OPTIONS_BUILD ${TEMP}
    CACHE STRING
    "A string used to set paths in platforms directories (not user-settable):
    \${HELYX_PRECISION_OPTION}Int\${HELYX_LABEL_SIZE}"
    FORCE
    )
mark_as_advanced(HELYX_OPTIONS_BUILD)


# ------------------------------  Output paths ------------------------------- #

# It is not possible to control the cmake build directory location from within
# CMake.  By the time these options are parsed, CMAKE_BINARY_DIR has already
# been set.  Re-setting this variable is advised against, and having tried it I
# can see why.  Changing the build directory paths must be done by the thing
# that calls CMake (e.g. emake, IDE's, etc...).

# There's a scenario where HELYX_OPTIONS changes and triggers a re-
# configuration.  Unless this is first caught by emake or an IDE, this re-
# configuration will take place in the wrong directory, and build progress will
# be wiped out for both the old and new configurations.
set(STRICT_BINARY_DIR_CHECKING
    ON
    CACHE BOOL
    "Strict check to ensure that CMAKE_BINARY_DIR matches HELYX_OPTIONS"
    FORCE
    )
mark_as_advanced(FORCE STRICT_BINARY_DIR_CHECKING)


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib
    CACHE PATH
    "Output directory for HELYX libraries"
    FORCE
    )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/bin
    CACHE PATH
    "Output directory for HELYX libraries"
    FORCE
    )
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/staticLib
    CACHE PATH
    "Output directory for static libraries (normally unused)"
    FORCE
    )
mark_as_advanced(FORCE
    CMAKE_LIBRARY_OUTPUT_DIRECTORY
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    )

# This variable left unforced, as it's often set on the command line (e.g. by
# CI systems)
set(CMAKE_INSTALL_PREFIX
    "${HELYX_PROJECT_DIR}/distribution/"
    CACHE STRING
    "The path used for the install target"
    )


# ============================================================================ #
# ----------------------------------  Misc ----------------------------------- #
# ============================================================================ #

option(FORCE_COLOURED_OUTPUT
    "Coloured output from compilers (do not use when piping to file)"
    OFF
)

option(FORCE_COLOURED_CMAKE_MESSAGES
    "Coloured output from CMake and MakeFiles (do not use when piping to file)"
    OFF
)

set(CMAKE_EXPORT_COMPILE_COMMANDS
    OFF
    CACHE BOOL
    "Enable/Disable output of compile commands during generation (used for clangd)"
    FORCE
    )

set(CMAKE_VERBOSE_MAKEFILE
    OFF
    CACHE BOOL
    "Debug build commands"
    FORCE
    )

set(HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL
    OFF
    CACHE BOOL
    "Switch to enable developers to set different optimisation levels for different targets"
    FORCE
    )
mark_as_advanced(FORCE HELYX_ENABLE_VARIABLE_OPTIMISATION_LEVEL)

set(HELYX_ADDITIONAL_MODULES
    ""
    CACHE STRING
    "Extra HELYX modules to look for in $HELYX_PROJECT_DIR/modules"
    FORCE
    )
mark_as_advanced(FORCE HELYX_ADDITIONAL_MODULES)

# The usage of lnIncludes cause minor (and not-so-minor) problems during compilation, for example:
#   • Generating files in the source dir is bad practice,
#   • There cannot be more than one file with the same name in the same target,
#   • The line “#include Foo.H” is ambiguous, and if Foo.H exists for many different targets then the
#     Foo.H that gets included can change depending on the order in which those targets are included in
#     the build command will fail to clean these directories properly, either with warnings or errors.
#     This is essentially unfixable by Engys without removing lnIncludes.
#
# For the reasons above, lnIncludes started to be deprecated in HELYX version 4.0.0
# and they will be completely deprecated on the next releases.
# HELYX-Core and Modules no longer relies on the lnIncludes folders, however, custom codes can still be compiled using the lnIncludes.
#
# LnIncludes will be completely deprecated on the subsequent releases.
# As one of the main changes, the helyx_include_directories() function (mostly used for custom code compilation
# against Helyx) will become obsolete and, at that point, users will need to use the helyx_additional_includes() instead.
#
# More details can be found in the HELYX-Core compilation Guide.
# To activate the generation of lnInclude folders, set HELYX_GENERATE_LNINCLUDES to ON.
set(HELYX_GENERATE_LNINCLUDES
    OFF
    CACHE BOOL
    "Generate lnInclude directories"
    FORCE
)
mark_as_advanced(FORCE HELYX_GENERATE_LNINCLUDES)

# Users are encouraged in the strongest possible terms to use bash or some
# other POSIX-derived shell.  Please only use csh if you have no other choice.
# Experimental support for python is also available.  Multiple different shells
# can be specified by setting RUNTIME_SHELL as a ;-separated list.  For example,
# to generate .shrc, .cshrc, and .py files you would specify "sh;csh;python".
set(RUNTIME_SHELL
    "sh"
    CACHE STRING
    "Shell for which to generate sourceable files (\"sh\", \"csh\", or \"python\")."
    FORCE
    )
set_property(CACHE RUNTIME_SHELL PROPERTY STRINGS sh csh python)
mark_as_advanced(FORCE RUNTIME_SHELL)

# These regexes are combined with the "|" specifier, signifying "or".  To
# disable all warnings, set ENVIRONMENT_VARIABLE_WARNING_WHITELIST to ".*".
# Note that these regexes must comply with the CMake regex specification, which
# is quite restrictive.
# At the time of writing, the CMake regex specification can be found here:
#   https://cmake.org/cmake/help/latest/command/string.html#regex-specification
set(ENVIRONMENT_VARIABLE_WARNING_WHITELIST
    "^FOAM_CONFIG;^FOAM_SIGFPE;^FOAM_USER_APPBIN;^FOAM_USER_LIBBIN"
    CACHE STRING
    "Semicolon-separated list of regexes to ignore when warning about environment variables"
    FORCE
    )
mark_as_advanced(FORCE ENVIRONMENT_VARIABLE_WARNING_WHITELIST)


# To also pack mpi-dependend library of an additional MPI build
# This is basically the HELYX_MPI_NAME of the additional build
# **/HELYXcore/platforms/<HELYX_OPTIONS>/lib/<HELYX_MPI_NAME>
set(HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS
    ""
    CACHE STRING
    "To pack additional MPI-dependent libraries from **/platforms/<HELYX_OPTIONS>/lib/<HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS>"
    )
mark_as_advanced(HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS)


# ============================================================================ #
# ---------------------  System specific compile flags  ---------------------- #
# ============================================================================ #

set(DEFERRED_DIRECTORY_DEFINITIONS
    ""
    CACHE STRING
    "A pair of directory and flags for specific compile definitions of that directory.
    The first entry in the directly were the options will be applied, and the second entry in a list of flags separated by space"
    FORCE
    )
mark_as_advanced(FORCE DEFERRED_DIRECTORY_DEFINITIONS)
set(DEFERRED_DIRECTORY_OPTIONS
    ""
    CACHE STRING
    "A pair of directory and flags for specific compile options of that directory.
    The first entry in the directly were the options will be applied, and the second entry in a list of flags separated by space"
    FORCE
    )
mark_as_advanced(FORCE DEFERRED_DIRECTORY_OPTIONS)


# ============================================================================ #
# --------------------------- Build Optimizations  --------------------------- #
# ============================================================================ #

# CMake Unity build is only available for CMake version >= 3.16
if ("${CMAKE_VERSION}" VERSION_LESS "3.16" OR "${CMAKE_CXX_COMPILER}" MATCHES mingw)
    set(CMAKE_UNITY_BUILD
        OFF
        CACHE BOOL
        "Enable unity builds"
        )
else ()
    set(CMAKE_UNITY_BUILD
        ON
        CACHE BOOL
        "Enable unity builds"
        )
endif()
mark_as_advanced(CMAKE_UNITY_BUILD)
set(CMAKE_UNITY_BUILD_BATCH_SIZE
    90
    CACHE STRING
    "Max number of files to unify in unity builds"
    FORCE
    )
mark_as_advanced(FORCE CMAKE_UNITY_BUILD_BATCH_SIZE)
set(TARGETS_TO_EXCLUDE_FROM_UNITY
    ""
    CACHE STRING
    "List of targets to exclude from the Unity build (useful for faster re-compilation when developing)"
    FORCE
    )
mark_as_advanced(FORCE TARGETS_TO_EXCLUDE_FROM_UNITY)

# There's literally no reason you'd ever want to turn this off, because depending
# on inline symbols is undefined anyway :D.
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# LTO defaults to ON in optimised builds only.
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Opt")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON CACHE BOOL "")
else()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF CACHE BOOL "")
endif()

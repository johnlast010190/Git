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
    (c) 2019 Engys Ltd.

Description
    CMake configuration options that don't directly affect the build.

[----------------------------------------------------------------------------]]

cmake_minimum_required(VERSION 3.13)

# Allow use of the IN_LIST operator
cmake_policy(SET CMP0057 NEW)

# This policy determines whether the list command will ignore empty
# elements in the list.  CMake 2.4 and below list commands ignored all
# empty elements in the list.  For example, a;b;;c would have length 3
# and not 4.  The OLD behavior for this policy is to ignore empty list
# elements.  The NEW behavior for this policy is to correctly count
# empty elements in a list.
# This policy was introduced in CMake version 2.6.0.  CMake version
# 3.10.2 warns when the policy is not set and uses OLD behavior.
cmake_policy(SET CMP0007 NEW)


set(HELYX_DIR_NAME
    "${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}"
    CACHE INTERNAL
    "The name of HELYX top-level directory"
    )

set(THIRDPARTY_DIR_NAME
    "ThirdParty-${HELYX_THIRDPARTY_VERSION}"
    CACHE INTERNAL
    "The name of ThirdParty top-level directory"
    )

# Add path for custom findModules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/findModules)


# Check generators version
if(CMAKE_GENERATOR MATCHES "Unix Makefiles" )
    # minimum suported Make version
    set(MAKE_VERSION_MIN "4.0.0")
    find_package(Make QUIET)
elseif(CMAKE_GENERATOR MATCHES "Ninja" )
    # test ninja and warn if ninja version is 1.10
    find_package(Ninja QUIET)
endif()

# CMake does a good job of automatically enabling/disabling coloured output based
# on whether it is connected to a terminal. Unfortunately, emake ruins that, so we
# have emake check if *it* is connected to an interactive terminal and tell us here.
# This allows us to default to using coloured output when users run emake in an
# interactive terminal, but automatically turn it off when emake output is redirected
# to a file or a pipe. The FORCE_* options continue to be an override.
if ("${CONNECTED_TO_TERMINAL}" OR "${FORCE_COLOURED_OUTPUT}")
    if (FORCE_COLOURED_OUTPUT)
        message(STATUS
            "Forcing coloured compiler output (if output is piped to file, colour information may also be written to file)"
        )
    endif()

    # Enable colours for compiler output
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()


if ("${CONNECTED_TO_TERMINAL}" OR "${FORCE_COLOURED_CMAKE_MESSAGES}")
    if (FORCE_COLOURED_CMAKE_MESSAGES)
        message(STATUS
            "Forcing coloured CMake messages (if output is piped to file, colour information may also be written to file)"
        )
    endif()

    # Enable coloured MakeFile (if used)
    if (CMAKE_GENERATOR MATCHES "(M|m?)ake$")
        message("CMAKE_COLOR_MAKEFILE=ON")
        set(CMAKE_COLOR_MAKEFILE ON)
    endif()
    # Colours definition
    if(NOT WIN32)
        string(ASCII 27 Esc)
        set(ColourReset       "${Esc}[m")
        set(ColourBold        "${Esc}[1m")
        set(ColourRed         "${Esc}[31m")
        set(ColourGreen       "${Esc}[32m")
        set(ColourYellow      "${Esc}[33m")
        set(ColourBlue        "${Esc}[34m")
        set(ColourMagenta     "${Esc}[35m")
        set(ColourCyan        "${Esc}[36m")
        set(ColourWhite       "${Esc}[37m")
        set(ColourBoldRed     "${Esc}[1;31m")
        set(ColourBoldGreen   "${Esc}[1;32m")
        set(ColourBoldYellow  "${Esc}[1;33m")
        set(ColourBoldBlue    "${Esc}[1;34m")
        set(ColourBoldMagenta "${Esc}[1;35m")
        set(ColourBoldCyan    "${Esc}[1;36m")
        set(ColourBoldWhite   "${Esc}[1;37m")
    endif()
else()
    # defaults for not coloured output
    set(CMAKE_COLOR_MAKEFILE OFF)
    if(NOT WIN32)
        string(ASCII 27 Esc)
        set(ColourReset       "")
        set(ColourBold        "")
        set(ColourRed         "")
        set(ColourGreen       "")
        set(ColourYellow      "")
        set(ColourBlue        "")
        set(ColourMagenta     "")
        set(ColourCyan        "")
        set(ColourWhite       "")
        set(ColourBoldRed     "")
        set(ColourBoldGreen   "")
        set(ColourBoldYellow  "")
        set(ColourBoldBlue    "")
        set(ColourBoldMagenta "")
        set(ColourBoldCyan    "")
        set(ColourBoldWhite   "")
    endif()
endif()

if (${CMAKE_EXPORT_COMPILE_COMMANDS})
#     message(STATUS
#         "Copying compile commands to project root for clangd integration"
#     )
    # At this point, compile_commands.json is not yet generated. Therefore, symlink
    execute_process(
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/compile_commands.json
        COMMAND ln -sf ${PROJECT_BINARY_DIR}/compile_commands.json ${PROJECT_SOURCE_DIR}/compile_commands.json
        TIMEOUT 10
        RESULT_VARIABLE copy_result
        OUTPUT_QUIET
    )
    if(NOT "0" STREQUAL "${copy_result}")
        message(WARNING "Failed to copy compile_commands.json")
    endif()
endif()
# Default context for files not in scope for clion
add_custom_target(.defaultContext)

list(APPEND CMAKE_LIBRARY_PATH
${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}
)
# This seems necessary to find the Helyx libraries for out-of-tree builds
# (not sure why CMAKE_LIBRARY_PATH above is not sufficient)
link_directories(
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}
)

# Set the RPATH to a relative path so that we don't have to rely on re-setting
# the LD_LIBRARY_PATH correctly, but keep library location flexible
# This is necessary for libraries that link to libraries - then the
# above CMAKE_LIBARARY_PATH does not help to locate them
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../lib/${HELYX_MPI_NAME}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The following settings are for packaging
# These settings must be set before adding targets, as target install logic
# happens when the target is added, to maintain backwards compatability with
# CMake < 3.13.
# Note that the COMPONENT name is appended to CPACK_PACKAGE_FILE_NAME to create
# the name of the package output.  Note that every component must have a name,
# which leads to this slightly strange naming convention.
# The package will be called ${CPACK_PACKAGE_FILE_NAME}-${COMPONENT_NAME}, and
# ${COMPONENT_NAME} cannot be empty.  There's probably an elegant way around
# this for Core, but for now I'll just call everything Windows-<component>,
# apart from core which just gets called "Windows".
if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    set(HELYX_PACKAGE_NAME "Windows"
        CACHE INTERNAL "Display name for packages (e.g. Ubuntu18, Windows)"
        FORCE
    )
elseif("${HELYX_PACKAGE_NAME}" STREQUAL "")
    set(HELYX_PACKAGE_NAME "POSIX" CACHE STRING
        "Display name for packages (e.g. Ubuntu18, Windows)"
        FORCE
    )
endif()
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME
    "${HELYX_PACKAGE_NAME}"
    CACHE INTERNAL
    "The default component name (POSIX or MSwindows)"
    )


# Here's a little widget to make generating graphviz graphs easier
add_custom_target(generateGraphsDir
    COMMAND ${CMAKE_COMMAND} -E echo "To customise the graphviz output, see instructions on the following page:"
    COMMAND ${CMAKE_COMMAND} -E echo "    https://cmake.org/cmake/help/latest/module/CMakeGraphVizOptions.html#variables-specific-to-the-graphviz-support"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${HELYX_PROJECT_DIR}/graphs"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${HELYX_PROJECT_DIR}/graphs/dotFiles"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${HELYX_PROJECT_DIR}/graphs/svgFiles"
)
add_custom_target(generateDotFiles
    COMMAND ${CMAKE_COMMAND} "--graphviz=HELYX.dot" ${CMAKE_BINARY_DIR}
    DEPENDS generateGraphsDir
    WORKING_DIRECTORY "${HELYX_PROJECT_DIR}/graphs/dotFiles"
)
add_custom_target(generateGraphs
    COMMAND dot -Tsvg ${HELYX_PROJECT_DIR}/graphs/dotFiles/HELYX.dot -o HELYX.svg
    DEPENDS generateDotFiles
    WORKING_DIRECTORY "${HELYX_PROJECT_DIR}/graphs/svgFiles"
)


# Strictly speaking, we should link against .dll.a files.  However, most of our
# thirdParty libraries only come with .dll files, and they work fine.  ".dll"
# has been removed from CMAKE_FIND_LIBRARY_SUFFIXES in very modern versions of
# CMake, so this just guarantees that our .dlls will be found.
if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
endif()


# We need to keep a list of unresolved links to include targets, in order to
# catch typos at configuration-time rather than link-time
set(UNRESOLVED_LINKS_TO_INCLUDES ""
    CACHE INTERNAL
    ""
)

#  Now tell people what's going on
message("")
message(TITLE "HELYX Configuration")
message(CLEAN "|     ${ColourBoldRed}  o  ${ColourReset}      |                                                             |
|  ${ColourBoldRed}  o     o  ${ColourReset}   |  HELYX (R) : Open-source CFD for Enterprise                 |")
message(PIPE_TERMINATED  "| ${ColourBoldRed}  o   O   o  ${ColourReset}  |  Version : ${HELYX_PROJECT_VERSION}")
message(CLEAN "|  ${ColourBoldRed}  o     o  ${ColourReset}   |  ENGYS Ltd. <http://engys.com/>                             |
|     ${ColourBoldRed}  o  ${ColourReset}      |                                                             |
[------------------------------------------------------------------------------]")
message(PIPE_TERMINATED  "| Build platform | ${HELYX_BUILD_PLATFORM}  ")
message(PIPE_TERMINATED  "| Compiler       | ${HELYX_COMPILER_NAME}")
message(PIPE_TERMINATED  "| Precision      | ${HELYX_PRECISION_OPTION}  ")
message(PIPE_TERMINATED  "| Label size     | Int${HELYX_LABEL_SIZE}")
message(PIPE_TERMINATED  "| Build type     | ${CMAKE_BUILD_TYPE}")
message(PIPE_TERMINATED  "| Target OS type | ${HELYX_SYSTEM_NAME}")
message(CLEAN "[------------------------------------------------------------------------------]")
message(PIPE_TERMINATED  "| CMake version  | ${CMAKE_VERSION}")
message(PIPE_TERMINATED  "| Generator      | ${CMAKE_GENERATOR}")
message(CLEAN
"[==============================================================================]

")


# Check the minimum supported version of GNU/Clang compiler
string(CONCAT compiler_message
"The minimum version of GNU/Clang compiler supported by HELYX is version 7. "
"However, the following compiler was detected: "
"${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}. "
"Please upgrade your GNU/Clang compiler to version 7+. \n"
)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 7.0.0)
        message(FATAL_ERROR ${compiler_message})
    endif()
endif()

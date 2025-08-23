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

[----------------------------------------------------------------------------]]

set(TARGET_NAME "OSspecific")

add_library(${TARGET_NAME} OBJECT
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigFpe.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigSegv.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigInt.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigQuit.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigStopAtWriteNow.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/signals/sigWriteNow.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/timer.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/fileStat.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/fileMonitor.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/regExp.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/cpuTime/cpuTime.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/clockTime/clockTime.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/memInfo/memInfo.C
${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/${HELYX_SYSTEM_NAME}.C
)

# Append to list of Helyx targets for emake autocompletion
get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
list(APPEND TARGETS ${TARGET_NAME})
set_property(GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS ${TARGETS})

if("POSIX" STREQUAL "${HELYX_SYSTEM_NAME}")
    if (DEFINED ENV{SunOS64})
        set(PRINTSTACK_SOURCE dummyPrintStack.C)
    else ()
        set(PRINTSTACK_SOURCE printStack.C)
    endif ()
    target_sources(${TARGET_NAME}
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/cpuInfo/cpuInfo.C
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/clockValue/clockValue.C
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/${PRINTSTACK_SOURCE}
    )

    find_file(BFD_HEADER NAMES "bfd.h")
    find_library(BFD_LIBRARY NAMES "libbfd.so" "libbfd.so.2")
    if (NOT EXISTS BFD_LIBRARY)
        # We don't want to use the static library, because it might have extra
        # unpredictable dependencies that we don't check for. However, use its
        # location as a hint to find the .so which (on OpenSUSE) has a super
        # complicated name with a version number and date
        find_library(BFD_STATIC_LIBRARY NAMES "libbfd.a")
        if (BFD_STATIC_LIBRARY)
            get_filename_component(BFD_DIRECTORY "${BFD_STATIC_LIBRARY}" DIRECTORY)
            file(GLOB BFD_SO_LIBS RELATIVE "${BFD_DIRECTORY}" "${BFD_DIRECTORY}/libbfd-2.*.so")
            if (BFD_SO_LIBS)
                list(GET BFD_SO_LIBS 0 BFD_SO_LIB)
                find_library(BFD_LIBRARY NAMES "${BFD_SO_LIB}")
            endif ()
        endif ()
    endif ()
    if (EXISTS ${BFD_HEADER} AND BFD_LIBRARY)
        message(STATUS "Found bfd library and header - enabling source/line output in stack trace")
        target_compile_definitions(${TARGET_NAME} PUBLIC -DUSE_BFD)
        # N.B.  These libraries are linked in libOpenFOAM, because only CMake >
        # 3.12 allows linking in object libraries
        # (Also, suspect these should be INTERFACE)
        # target_link_libraries(${TARGET_NAME} ${BFD_LIBRARY})
    else ()
        message(STATUS "<bfd.h> or bfd library not found; disabling source/line output in stack trace")
    endif()

elseif("MSwindows" STREQUAL "${HELYX_SYSTEM_NAME}")
    target_sources(${TARGET_NAME}
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/printStack.C
    )
    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_REGEX}
        PRIVATE ${THIRDPARTY_STACK_TRACE}
    )
endif()

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/${HELYX_SYSTEM_NAME}/memInfo/memInfo.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )

add_helyx_include_target(${CMAKE_CURRENT_SOURCE_DIR}/${HELYX_SYSTEM_NAME})
add_to_dir_target(OSspecific)

find_file(INOTIFY_HEADER NAMES "sys/inotify.h")
if (EXISTS ${INOTIFY_HEADER})
    message(STATUS "Found <sys/inotify.h> - enabling inotify for file monitoring")
    target_compile_definitions(${TARGET_NAME} PUBLIC -DUSE_INOTIFY)
endif()

# In older versions of CMake, OBJECT libraries can't link to anything -> have to
# manually add cwd as include dir
target_include_directories(${TARGET_NAME}
    PRIVATE ${HELYX_PROJECT_DIR}/src/OpenFOAM
    PRIVATE ${HELYX_PROJECT_DIR}/src/Pstream/mpi
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${HELYX_SYSTEM_NAME}

    PRIVATE ${THIRDPARTY_MPI_INC}
)

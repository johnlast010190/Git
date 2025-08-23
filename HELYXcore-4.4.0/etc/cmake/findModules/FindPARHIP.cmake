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
    (c) 2021 Engys Ltd.

Description
Custom findModule for ParHIP

[----------------------------------------------------------------------------]]

set(PARHIP_PATH_SUFFIXES "../lib/${HELYX_MPI_NAME};../lib;../${HELYX_MPI_NAME};${PARHIP_PATH_SUFFIXES}")
helyx_find_thirdparty_package(PARHIP "parhip_interface" "parhip_interface.h")

if(PARHIP_FOUND)
    find_program (PARHIP_EXE
        NAMES parhip
        HINTS "${PARHIP_ARCH_PATH}"
        PATH_SUFFIXES "bin" "${PARHIP_PATH_SUFFIXES}"
        NO_CACHE
        )

    string(CONCAT determinism_warning
        "It was therefore not possible to detect whether ParHIP was "
        "compiled with the DETERMINISTIC_PARHIP flag.  Note that this "
        "flag is only available in ParHIP version 3.14 and later.\n"
        "Use of non-determistic code may make "
        "it difficult to get reproducible results from chaotic cases.\n"
        "HELYXcore will still attempt to compile with ParHIP support."
    )
    if(PARHIP_EXE)
        execute_process(COMMAND ${PARHIP_EXE} --version
            TIMEOUT 1
            RESULT_VARIABLE return_code
            ERROR_VARIABLE parhip_error_string
            OUTPUT_VARIABLE parhip_version_string
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if (NOT "${return_code}" STREQUAL 0)
            string(CONCAT s
                "The following command returned a non-zero error code:\n"
                "\t${PARHIP_EXE} --version\n"
                "${determinism_warning}"
            )
            message(WARNING "${s}")
        elseif (NOT "${parhip_error_string}" STREQUAL "")
            string(CONCAT s
                "\"parhip --version\" printed the following error:\n"
                "\t${parhip_error_string}\n"
                "${determinism_warning}"
            )
            message(WARNING "${s}")
        elseif(NOT "${parhip_version_string}" MATCHES "DETERMINISTIC_PARHIP")
            string(CONCAT s
                "ParHIP reported the following version information:\n"
                "\t${parhip_version_string}\n"
                "This string does not contain the substring "
                "\"DETERMINISTIC_PARHIP\".  It therefore appears that this "
                "version of ParHIP was not compiled with the "
                "DETERMINISTIC_PARHIP flag (available since ParHIP version "
                "3.14)."
                "Use of non-determistic code may make "
                "it difficult to get reproducible results from chaotic cases.\n"
                "HELYXcore will still attempt to compile with ParHIP support."
            )
            message(WARNING "${s}")
        endif()
    else()
        string(CONCAT s
            "The ParHIP library was found, but the \"parhip\" executable was "
            "not."
            "${determinism_warning}"
        )
        message(WARNING "${s}")
    endif()
endif()

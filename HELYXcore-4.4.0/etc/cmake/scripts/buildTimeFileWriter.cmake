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
    A cmake script for writing simple files at build time.  It must be called
    with two arguments, "filename" and "content", as shown:
        cmake -Dfilename=test_file.txt -Dcontent='Hello world!'  \
        -P ${CMAKE_SOURCE_DIR}/etc/cmake/scripts/buildTimeFileWriter.cmake

[----------------------------------------------------------------------------]]


cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

foreach(expected_variable filename content)
    if(NOT ${expected_variable})
        message(SEND_ERROR "BuildTimeFileWriter.cmake called with incorrect parameters - ${expected_variable} evaluated to false")
    endif()
endforeach()

file(WRITE ${filename} ${content})

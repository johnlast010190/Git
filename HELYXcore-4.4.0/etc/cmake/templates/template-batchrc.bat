@echo off
rem [--------------------------------------------------------------------------]
rem |       o        |
rem |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
rem |   o   O   o    |  Version : @HELYX_PROJECT_VERSION@
rem |    o     o     |  ENGYS Ltd. <http://engys.com/>
rem |       o        |
rem [--------------------------------------------------------------------------]
rem License
rem
rem    This file is part of HELYXcore.
rem    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
rem
rem    HELYXcore is free software: you can redistribute it and/or modify it
rem    under the terms of the GNU General Public License as published by
rem    the Free Software Foundation, either version 3 of the License, or
rem    (at your option) any later version.
rem
rem    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
rem    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
rem    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
rem    for more details.
rem
rem    You should have received a copy of the GNU General Public License
rem    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.
rem
rem Copyright
rem    (c) 2011-2013  blueCape Lda.
rem    (c) 2019-2020 Engys Ltd.
rem
rem  Script
rem      batchrc.bat
rem
rem  Description
rem    This batch file initializes the OpenFOAM environment for running in the
rem    Windows Command Line.
rem
rem [--------------------------------------------------------------------------]

rem Isolate HELYX from other applications
set PATH=%WINDIR%;%WINDIR%\SYSTEM32

set CURR_DIR=%cd%

%~d0
cd "%~dp0.."
call :SETHOME "%CD%"

rem BEGINING of Default options --------------------------------------

rem - Compiler:
rem PICK ONE from these:
rem     mingw32    - custom build of the Gcc+mingw cross-compiler
rem     mingw_w32  - custom build of the 32bit Gcc+mingw-w64 cross-compiler
rem     mingw_w64  - custom build of the 64bit Gcc+mingw-w64 cross-compiler
set HELYX_COMPILER=mingw_w64

rem - MPI implementation:
rem PICK ONE: "", OPENMPI, MPICH, MSMPI2008, MSMPI2012
set HELYX_MPLIB=MSMPI80

rem - Precision:
rem PICK ONE: SP or DP
set HELYX_PRECISION_OPTION=DP

rem - label size
set HELYX_LABEL_SIZE=32

rem END of Default options --------------------------------------------

rem BEGGINING of Process summoned options ---------------------------
FOR %%A IN (%*) DO (
  set %%A
)
rem END of Process summoned options ---------------------------------

set USER=hxuser
set USERNAME=ofuser
set HELYX_PROJECT_NAME=HELYXcore
set HELYX_PROJECT_VERSION=@HELYX_PROJECT_VERSION@
set HELYX_THIRDPARTY_VERSION=@HELYX_THIRDPARTY_VERSION@
set HELYX_PROJECT_DIR=%HOME%\%HELYX_PROJECT_NAME%-%HELYX_PROJECT_VERSION%
set HELYX_PROJECT_NAME_USER_DIR=%HOME%\%USER%-%HELYX_PROJECT_VERSION%
set HELYX_THIRDPARTY_DIR=%HOME%\ThirdParty-%HELYX_THIRDPARTY_VERSION%

set HELYX_OS=MSwindows
set HELYX_BUILD_PLATFORM=linux
set HELYX_BUILD_TYPE=@CMAKE_BUILD_TYPE@
set HELYX_COMPILER_ARCH=""
set HELYX_COMPILER_LIB_ARCH=""

IF "%HELYX_COMPILER%"=="mingw32"   set HELYX_BUILD_PLATFORM_OPTION=32
IF "%HELYX_COMPILER%"=="mingw-w32" set HELYX_BUILD_PLATFORM_OPTION=32
IF "%HELYX_COMPILER%"=="mingw-w64" set HELYX_BUILD_PLATFORM_OPTION=64

IF "%HELYX_COMPILER%"=="mingw32"   set HELYX_COMPILER_ARCH=i686-pc-mingw32
IF "%HELYX_COMPILER%"=="mingw-w32" set HELYX_COMPILER_ARCH=i686-pc-mingw32
IF "%HELYX_COMPILER%"=="mingw-w64" set HELYX_COMPILER_ARCH=x86_64-w64-mingw32

rem - Floating-point signal handling:
rem     set or unset
set FOAM_SIGFPE=""

rem - memory initialisation:
rem     set or unset
rem set FOAM_SETNAN=""

set HELYX_C_COMPILER=%HELYX_COMPILER_ARCH%-gcc
set HELYX_CXX_COMPILER=%HELYX_COMPILER_ARCH%-g++
set FOAM_JOB_DIR=%HELYX_PROJECT_DIR%\..\jobControl

rem wmake configuration
set HELYX_DIR=%HELYX_PROJECT_DIR%\wmake
set HELYX_LINK_LANGUAGE=c++
set HELYX_OPTIONS=%HELYX_BUILD_PLATFORM%%HELYX_COMPILER%%HELYX_PRECISION_OPTION%%HELYX_BUILD_TYPE%

rem base executables/libraries
set HELYX_RUNTIME_OUTPUT_DIRECTORY=%HELYX_PROJECT_DIR%\platforms\%HELYX_BUILD_PLATFORM%%HELYX_COMPILER%Gcc%HELYX_PRECISION_OPTION%Int%HELYX_LABEL_SIZE%%HELYX_BUILD_TYPE%\bin
set HELYX_LIBRARY_OUTPUT_DIRECTORY=%HELYX_PROJECT_DIR%\platforms\%HELYX_BUILD_PLATFORM%%HELYX_COMPILER%Gcc%HELYX_PRECISION_OPTION%Int%HELYX_LABEL_SIZE%%HELYX_BUILD_TYPE%\lib

rem external (ThirdParty) libraries
set FOAM_EXT_LIBBIN=%HELYX_THIRDPARTY_DIR%\platforms\%HELYX_BUILD_PLATFORM%%HELYX_COMPILER%Gcc%HELYX_PRECISION_OPTION%Int%HELYX_LABEL_SIZE%%HELYX_BUILD_TYPE%\lib

rem user executables/libraries
set HELYX_USER_APPBIN=%HELYX_PROJECT_USER_DIR%\platforms\%HELYX_OPTIONS%\bin
set HELYX_USER_LIBBIN=%HELYX_PROJECT_USER_DIR%\platforms\%HELYX_OPTIONS%\lib

IF "%HELYX_MPLIB%"=="""" set HELYX_MPI_NAME=dummy
IF "%HELYX_MPLIB%"=="OPENMPI" set HELYX_MPI_NAME=openmpi-1.6.2
IF "%HELYX_MPLIB%"=="MSMPI2008" set HELYX_MPI_NAME=msmpi-2008R2
IF "%HELYX_MPLIB%"=="MSMPI71" set HELYX_MPI_NAME=MS-MPI-7.1
IF "%HELYX_MPLIB%"=="MSMPI80" set HELYX_MPI_NAME=MS-MPI-8.0

set MPI_HOME=%HELYX_THIRDPARTY_DIR%\%HELYX_MPI_NAME%
set MPI_ARCH_PATH=%HELYX_THIRDPARTY_DIR%\platforms\%HELYX_BUILD_PLATFORM%%HELYX_COMPILER%Gcc%HELYX_PRECISION_OPTION%Int%HELYX_LABEL_SIZE%\%HELYX_MPI_NAME%

IF "%HELYX_MPLIB%"=="OPENMPI" set OPAL_PKGDATADIR=%MPI_ARCH_PATH%\share\openmpi
IF "%HELYX_MPLIB%"=="MPICH" set MPICH_ROOT=%MPI_ARCH_PATH%

set HELYX_MPI_NAME_LIBBIN=%HELYX_LIBRARY_OUTPUT_DIRECTORY%\%HELYX_MPI_NAME%
set MPI_BUFFER_SIZE=20000000
set HX_MPI_ARCH_PATH="C:\Program Files\Microsoft MPI"
set HELYX_TUTORIALS=%HELYX_PROJECT_DIR%\examples

set PATH=%PATH%;%HELYX_DIR%;%HX_MPI_ARCH_PATH%\lib;%HX_MPI_ARCH_PATH%\bin;%GSL_ARCH_PATH%\lib;%HELYX_MPI_NAME_LIBBIN%;%HELYX_RUNTIME_OUTPUT_DIRECTORY%;%HELYX_LIBRARY_OUTPUT_DIRECTORY%;%FOAM_EXT_LIBBIN%;%HELYX_PROJECT_DIR%\bin;%HELYX_PROJECT_DIR%\..\..\GUI\ext\EVTK-OpenGL2;%HELYX_PROJECT_DIR%\..\..\GUI\ext\OCCT\lib;

rem add (non-dummy) MPI implementation
rem dummy MPI already added to LD_LIBRARY_PATH and has no external libraries
IF NOT "%HELYX_MPI_NAME%"=="dummy" set PATH=%PATH%;%FOAM_EXT_LIBBIN%\%HELYX_MPI_NAME%

rem Source all *.bat files present at "%HELYX_PROJECT_DIR%\etc\config.d"
for %%A in ("%HELYX_PROJECT_DIR%\etc\config.d\*.bat") DO CALL "%%A"
GOTO END

:SETHOME
set HOME=%~dp1
set HOME=%HOME:~0,-1%
GOTO :EOF

:END
rem cd "%HELYX_PROJECT_DIR/..%"
cd /D %CURR_DIR%

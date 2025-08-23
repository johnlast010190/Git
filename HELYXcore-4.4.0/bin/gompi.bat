@echo off
rem ------------------------------------------------------------------------------
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
rem    (c) 2011-2014  blueCape Lda.
rem
rem  Script
rem      gompi.bat
rem 
rem  Description
rem      This batch file is analogous to foamJob, but for Windows Command Line.
rem 
rem ------------------------------------------------------------------------------

set /A x_numprocs=0
set MACHINEFILE=
set PREVIOUSPATH=%CD%

rem Detect if the execution folder should be changed
set PARAMS= 
set CASEPATH= 
:LOOP
  if "%~1"=="" goto afterloop
  if "%~1"=="-case" set CASEPATH=%2
  set PARAMS=%PARAMS% %1
  shift
  goto LOOP
:AFTERLOOP

rem go into case folder and launch application with the proper MPI setup
call :GOTODRIVE %CASEPATH%
cd %CASEPATH%
call :CONTINUE %PARAMS%
cd %PREVIOUSPATH%
GOTO END

:GOTODRIVE 
%~d1
GOTO END

:CONTINUE
rem simple count
rem for /D %%a in (processor*) do @set /A x_numprocs=x_numprocs+1

rem read directly from decomposeParDict and break on first find
setlocal enableextensions enabledelayedexpansion
FOR /F "eol=» tokens=1,2 delims=; " %%i in (system\decomposeParDict) do (
  IF "%%i" == "numberOfSubdomains" (
    set x_numprocs=%%j
    GOTO :BREAKLOOP
    )
)
:BREAKLOOP

IF EXIST "hostfile" set MACHINEFILE=hostfile
IF EXIST "machines" set MACHINEFILE=machines
IF EXIST "system\hostfile" set MACHINEFILE=system\hostfile
IF EXIST "system\machines" set MACHINEFILE=system\machines

rem Generate temporary batch file name, which is necessary in cases with long path names
set RNDBATCH=%TIME:~6,2%%RANDOM%.bat

IF NOT "%MACHINEFILE%" == "" set MACHINEFILE=-machinefile %MACHINEFILE%

rem Due to Windows user account control restrictions, executables cannot
rem contain the word "patch".  These executables are renamed at compile 
rem time, and we must compensate for that here.
rem N.B. This is case-independent.
set COMMAND_NAME=%1
set COMMAND_NAME=%COMMAND_NAME:patch=ptch%

@echo on
mpiexec -n %x_numprocs% %MPI_ACCESSORY_OPTIONS% %MACHINEFILE% %COMMAND_NAME% -parallel %2 %3 %4 %5 %6 %7 %8 %9
@echo off
GOTO END


:END
set x_numprocs=
IF EXIST "%RNDBATCH%" del %RNDBATCH%

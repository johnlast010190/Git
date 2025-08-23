/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
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
    (c) 2017 OpenCFD Ltd.
    (c) 2011 OpenFOAM Foundation

Namespace
    Foam

Description
    Namespace for OpenFOAM


Global
    Foam::FOAMversion

Description
    OpenFOAM version number static string.

Global
    Foam::FOAMbuild

Description
    OpenFOAM build information as a static string

Global
    Foam::FOAMbuildArch

Description
    OpenFOAM build architecture information (endian, label/scalar sizes)
    as a static string

SourceFiles
    global.Cver

\*---------------------------------------------------------------------------*/

#ifndef foamVersion_H
#define foamVersion_H

#include <string>

// The directory name for user-resources (located in the HOME directory)
#define HELYX_USER_RESOURCE_DIRNAME  ".HELYXcore"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    extern const char* const FOAMversion;
    extern const char* const FOAMbuild;
    extern const std::string FOAMbuildArch;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

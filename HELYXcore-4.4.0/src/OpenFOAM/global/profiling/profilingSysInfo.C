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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "global/profiling/profilingSysInfo.H"
#include "global/foamVersion.H"
#include "global/clock/clock.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "include/OSspecific.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope function
inline static void printEnv
(
    Foam::Ostream& os,
    const Foam::word& key,
    const std::string& envName
)
{
    const std::string value = Foam::getEnv(envName);
    if (!value.empty())
    {
        os.writeEntry(key, value);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingSysInfo::profilingSysInfo()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingSysInfo::~profilingSysInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::profilingSysInfo::write
(
    Ostream& os
) const
{
    os.writeEntry("host",       Foam::hostName(false)); // short name
    os.writeEntry("date",       Foam::clock::dateTime());

    // compile-time information
    os.writeEntry("version",    std::string(FOAMversion));
    os.writeEntry("build",      std::string(FOAMbuild));

    printEnv(os, "arch",         "HELYX_BUILD_PLATFORM");
    printEnv(os, "compiler",     "HELYX_COMPILER");
//     printEnv(os, "mplib",        "WM_MPLIB");
    printEnv(os, "options",      "HELYX_OPTIONS");

    return os;
}


// ************************************************************************* //

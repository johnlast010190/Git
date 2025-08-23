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
    (c) 2016 OpenCFD Ltd.

Application

Description

\*---------------------------------------------------------------------------*/

#include "global/profiling/profilingSysInfo.H"
#include "db/IOstreams/IOstreams.H"
#include "primitives/endian/endian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    profiling::sysInfo().write(Info);

#ifdef HELYX_BIG_ENDIAN
    Info
        << "HELYX_BIG_ENDIAN is defined"
        << nl;
#endif
#ifdef HELYX_LITTLE_ENDIAN
    Info
        << "HELYX_LITTLE_ENDIAN is defined"
        << nl;
#endif

    Info<< "Runtime endian check: big=" << endian::isBig()
        << " little=" << endian::isLittle()
        << nl;

    return 0;
}


// ************************************************************************* //

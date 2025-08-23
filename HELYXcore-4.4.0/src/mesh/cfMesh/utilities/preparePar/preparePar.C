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
    (c) Creative Fields, Ltd.
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Application
    Prepares the case for a parallel mesh generation run

Description
    - creates processor* directories which contain data for processors

\*---------------------------------------------------------------------------*/

#include "db/objectRegistry/objectRegistry.H"
#include "global/argList/argList.H"
#include "db/Time/Time.H"

#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#include "include/setRootCase.H"
#include "include/createTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary decomposeParDict
    (
        IOobject
        (
            "decomposeParDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const label nProcessors
    (
        decomposeParDict.lookup<label>("numberOfSubdomains")
    );

    for (label procI=0;procI<nProcessors;++procI)
    {
        fileName file("processor");
        std::ostringstream ss;
        ss << procI;
        file += ss.str();
        Info<< "Creating " << file << endl;

        //- create a directory for processor data
        mkDir(runTime.path()/file);

        //- copy the contents of the const directory into processor*
        cp(runTime.path()/"constant", runTime.path()/file);

        //- generate 0 directories for
        mkDir(runTime.path()/file/"0");
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

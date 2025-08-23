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
    (c) 2011-2015 OpenFOAM Foundation

Application
    foamMeshToFluent

Group
    grpMeshConversionUtilities

Description
    Writes out the OpenFOAM mesh in Fluent mesh format.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fluentFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "include/setRootCase.H"
    #include "include/createTime.H"

    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    fluentFvMesh mesh
    (
        IOobject
        (
            fluentFvMesh::defaultRegion,
            runTime.constant(),
            runTime
        )
    );

    mesh.writeFluentMesh();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

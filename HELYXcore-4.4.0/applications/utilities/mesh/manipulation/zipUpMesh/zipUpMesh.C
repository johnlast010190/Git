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
    zipUpMesh

Group
    grpMeshManipulationUtilities

Description
    Reads in a mesh with hanging vertices and zips up the cells to guarantee
    that all polyhedral cells of valid shape are closed.

    Meshes with hanging vertices may occur as a result of split
    hex mesh conversion or integration or coupled math face pairs.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "polyMeshZipUpCells/polyMeshZipUpCells.H"
#include "primitives/bools/lists/boolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addRegionOption.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createNamedPolyMesh.H"

    if (polyMeshZipUpCells(mesh))
    {
        Info<< "Writing zipped-up polyMesh to " << mesh.facesInstance() << endl;
        mesh.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

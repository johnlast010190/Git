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
    attachMesh

Group
    grpMeshManipulationUtilities

Description
    Attach topologically detached mesh using prescribed mesh modifiers.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "polyTopoChange/attachPolyTopoChanger/attachPolyTopoChanger.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addOverwriteOption.H"
    argList::noParallel();

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    runTime.functionObjects().off();
    #include "include/createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");

    if (!overwrite)
    {
        runTime++;
    }

    Info<< "Time = " << runTime.timeName() << nl
        << "Attaching sliding interface" << endl;

    attachPolyTopoChanger(mesh).attach();

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Application
    orientFaceZone

Group
    grpMeshManipulationUtilities

Description
    Corrects orientation of faceZone.

    - correct in parallel - excludes coupled faceZones from walk
    - correct for non-manifold faceZones - restarts walk

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "orientFaceZone/orientFaceZone.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addRegionOption.H"
    argList::validArgs.append("faceZone");
    argList::validArgs.append("outsidePoint");

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createNamedPolyMesh.H"

    const word zoneName  = args[1];
    const point outsidePoint = args.argRead<point>(2);

    orientFaceZone orFaceZone(mesh, zoneName);
    orFaceZone.orientByPoint(outsidePoint);

    return 0;
}


// ************************************************************************* //

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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFormats/ac3d/AC3DsurfaceFormat.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// read MeshedSurface
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    AC3DsurfaceFormat,
    face,
    fileExtension,
    ac
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    AC3DsurfaceFormat,
    triFace,
    fileExtension,
    ac
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    AC3DsurfaceFormat,
    labelledTri,
    fileExtension,
    ac
);

// write MeshedSurfaceProxy
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    AC3DsurfaceFormat,
    face,
    write,
    fileExtension,
    ac
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    AC3DsurfaceFormat,
    triFace,
    write,
    fileExtension,
    ac
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    AC3DsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    ac
);


// write UnsortedMeshedSurface
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    AC3DsurfaceFormat,
    face,
    write,
    fileExtension,
    ac
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    AC3DsurfaceFormat,
    triFace,
    write,
    fileExtension,
    ac
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    AC3DsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    ac
);

}
}

// ************************************************************************* //

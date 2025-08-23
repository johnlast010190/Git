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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFormats/vtk/VTKsurfaceFormat.H"

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
    VTKsurfaceFormat,
    face,
    fileExtension,
    vtk
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    VTKsurfaceFormat,
    triFace,
    fileExtension,
    vtk
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    VTKsurfaceFormat,
    labelledTri,
    fileExtension,
    vtk
);

// write MeshedSurfaceProxy
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    VTKsurfaceFormat,
    face,
    write,
    fileExtension,
    vtk
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    VTKsurfaceFormat,
    triFace,
    write,
    fileExtension,
    vtk
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    VTKsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    vtk
);

// write UnsortedMeshedSurface
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    VTKsurfaceFormat,
    face,
    write,
    fileExtension,
    vtk
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    VTKsurfaceFormat,
    triFace,
    write,
    fileExtension,
    vtk
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    VTKsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    vtk
);

}
}

// ************************************************************************* //

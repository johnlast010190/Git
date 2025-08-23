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
    (c) 2012-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "surfaceFormats/vtk/VTKsurfaceFormat.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readVTK(const fileName& fName)
{
    // Read (and triangulate) point, faces, zone info
    fileFormats::VTKsurfaceFormat<triFace> surf(fName);

    List<labelledTri> tris(surf.size());
    forAll(tris, i)
    {
        const triFace& f = surf[i];
        tris[i] = labelledTri(f[0], f[1], f[2], 0);
    }

    // Add regions from zone
    const List<surfZone>& surfZones = surf.surfZones();

    geometricSurfacePatchList patches;

    if (surfZones.size())
    {
        patches.setSize(surfZones.size());
        forAll(surfZones, zoneI)
        {
            const surfZone& zone = surfZones[zoneI];

            // Add patch. Convert synthetic 'zone' name into 'patch' for now.
            // (vtk format does not contain region names)
            word regionName = zone.name();
            if (regionName != (string("zone") + name(zoneI)))
            {
                regionName = string("patch") + name(zoneI);
            }

            patches[zoneI] = geometricSurfacePatch
            (
                regionName,
                zoneI,
                zone.geometricType()
            );

            // Set triangle regions
            for (label i = zone.start(); i < zone.start()+zone.size(); ++i)
            {
                tris[i].region() = zoneI;
            }
        }
    }
    else
    {
        // Add single (default) patch
        // Triangle regions already set to 0
        patches = { geometricSurfacePatch("patch0", 0) };
    }


    // Create triSurface
    *this = triSurface
    (
        std::move(tris),
        patches,
        surf.stdMovePoints()
    );

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

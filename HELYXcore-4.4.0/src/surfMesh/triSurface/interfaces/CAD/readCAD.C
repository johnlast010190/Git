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
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "cad/CADReader.H"
#include "triSurface/triSurface.H"
#include "include/OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface
(
    const fileName& CADfileName,
    const scalar& linDeflection,
    const scalar& angDeflection,
    const bool forceTriangulation, /*= false*/
    const word nameRegionsBy, /*= "solid"*/
    const wordHashSet& ignoreRegions /*= wordHashSet()*/
)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr),
    points0Ptr_(nullptr)
{
    fileName triFile =
        CADfileName.lessExt()
      + "_OCCT_" + name(linDeflection)
      + "_" + name(angDeflection)
      + ".obj";

    if (isFile(triFile) && !forceTriangulation)
    {
        Info<< "Skipping triangulation."
            << nl << "Found " << triFile
            << endl;

        *this = triSurface(triFile);
    }
    else
    {
        fileFormats::CADReader reader
        (
            CADfileName,
            nameRegionsBy,
            /*distributed*/false,
            ignoreRegions
        );

        HashTable<label> patchIDs;

        reader.triangulate
        (
            linDeflection,
            angDeflection,
            storedPoints(),
            storedFaces(),
            patchIDs
        );

        patches_.setSize(patchIDs.size());

        forAllConstIters(patchIDs, iter)
        {
            const label patchIdx = iter.object();

            patches_[patchIdx] =
                geometricSurfacePatch
                (
                    iter.key(),
                    patchIdx
                );
        }

        if (Pstream::master())
        {
            Info<< "Writing OCCT triangulation to "
                << triFile
                << endl;
            write(triFile);
        }
    }
}

// ************************************************************************* //

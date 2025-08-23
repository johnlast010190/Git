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

Class
    triSurfaceRemoveFacets

Description
    Divides the surface mesh into regions bounded by feature edges

SourceFiles
    triSurfaceRemoveFacets.C
    triSurfaceRemoveFacetsFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceRemoveFacets_H
#define triSurfaceRemoveFacets_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceRemoveFacets Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceRemoveFacets
{
    // Private data

        //- Reference to triSurf
        triSurf& surf_;

        //- Patches/subsets for removal
        DynList<word> selectedEntities_;


    // Private member functions

        //- Remove facets in selected patches/subsets
        void markFacetsForRemoval(boolList&) const;


public:

    // Constructors

        //- Construct from octree
        triSurfaceRemoveFacets(triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfaceRemoveFacets(const triSurfaceRemoveFacets&) = delete;


    //- Destructor
    ~triSurfaceRemoveFacets();


    // Member Functions

        //- Add patch for removal
        void selectFacetsInPatch(const word&);

        //- Add subsets for removal
        void selectFacetsInSubset(const word&);

        //- Perform removal of selected facets
        void removeFacets();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceRemoveFacets&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

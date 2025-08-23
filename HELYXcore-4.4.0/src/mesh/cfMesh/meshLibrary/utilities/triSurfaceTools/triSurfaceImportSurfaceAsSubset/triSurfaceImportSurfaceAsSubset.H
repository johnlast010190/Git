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
    triSurfaceImportSurfaceAsSubset

Description
    Creates a subset in the master surface consisting of facets which are
    near the other surface

SourceFiles
    triSurfaceImportSurfaceAsSubset.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceImportSurfaceAsSubset_H
#define triSurfaceImportSurfaceAsSubset_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshOctree;

/*---------------------------------------------------------------------------*\
             Class triSurfaceImportSurfaceAsSubset Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceImportSurfaceAsSubset
{
    // Private data

        //- Reference to triSurf
        triSurf& surf_;

        //- Pointer to meshOctree, needed for searching on the master surface
        meshOctree* octreePtr_;


    // Private member functions

        void createOctree(const triSurf&, meshOctree&);


public:

    // Constructors

        //- Construct from triSurf
        triSurfaceImportSurfaceAsSubset(triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfaceImportSurfaceAsSubset
        (
            const triSurfaceImportSurfaceAsSubset&
        ) = delete;


    //- Destructor
    ~triSurfaceImportSurfaceAsSubset();


    // Member Functions

        //- finds the nearest faces in the surface to the import surf
        //- and creates a subset
        void addSurfaceAsSubset
        (
            const triSurf& importSurf,
            const word& subsetName,
            const scalar angleTol = 5.*M_PI/180.
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceImportSurfaceAsSubset&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

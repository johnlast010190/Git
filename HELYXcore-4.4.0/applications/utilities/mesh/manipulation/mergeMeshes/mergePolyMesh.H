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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

Class
    Foam::mergePolyMesh

Description
    Add a given mesh to the original mesh to create a single new mesh

SourceFiles
    mergePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mergePolyMesh_H
#define mergePolyMesh_H

#include "meshes/polyMesh/polyMesh.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class mergePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mergePolyMesh
:
    public polyMesh
{
    // Private data

        //- Topological change to accumulated all mesh changes
        polyTopoChange meshMod_;

        //- Patch names
        DynamicList<word> patchNames_;

        //- Patch dictionaries
        DynamicList<dictionary> patchDicts_;

        //- Point zone names
        DynamicList<word> pointZoneNames_;

        //- Face zone names
        DynamicList<word> faceZoneNames_;

        //- Cell zone names
        DynamicList<word> cellZoneNames_;


    // Private Member Functions

        //- Set patch index given a name and type
        void setPatchIndices
        (
            const polyBoundaryMesh& bm,
            labelList& patchIndices,
            DynamicList<labelPair>& movedPatches
        );

        //- Return zone index given a list of active zones and a name
        label zoneIndex(DynamicList<word>&, const word&);

        //- Shuffle processor patches to be last
        void sortProcessorPatches();


public:

    //- Runtime type information
    TypeName("mergePolyMesh");


    // Constructors

        //- Construct from IOobject
        mergePolyMesh(const IOobject& io);

        //- Construct from components (zero-sized mesh)
        mergePolyMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Disallow default bitwise copy construction
        mergePolyMesh(const mergePolyMesh&) = delete;

    //- Destructor
    virtual ~mergePolyMesh()
    {}


    // Member Functions

        //- Reset meshMod and whether to add to existing mesh
        void topoChangeMod(bool addMesh=true);

        //- Add a mesh
        void addMesh(const polyMesh& m);

        //- Merge meshes
        void merge();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mergePolyMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

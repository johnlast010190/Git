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
    (c) 2010 Engys Ltd.

Application
    Foam::addHexInternalLayers

Description
    Add layers of cells

\*---------------------------------------------------------------------------*/

#ifndef addHexInternalLayers_H
#define addHexInternalLayers_H

#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class polyTopoChangeMap;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class addHexInternalLayers Declaration
\*---------------------------------------------------------------------------*/

class addHexInternalLayers
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        meshDistribute& distributor_;

        List<labelPair> baffles_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("addHexInternalLayers");

    // Constructors

        //- Construct from mesh.
        addHexInternalLayers
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            meshDistribute& distributor
        );

    // Member Functions

            void addCells
            (
                const indirectPrimitivePatch& pp,
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& faceType,
                List<labelPair>& newFaceCells,
                List<labelPair>& newEdgeCells,
                List<labelPair>& newPointCells,
                polyTopoChange& meshMod
            );

            void addPoints
            (
                const vectorField& meshNormals,
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& baffleLevel,
                labelList& addedPoints,
                List<DynamicList<labelPair>>& convexEdgePts,
                List<DynamicList<labelPair>>& convexFacePts,
                polyTopoChange& meshMod
            );

            void calculateTypes
            (
                const indirectPrimitivePatch& pp,
                const vectorField& meshEdgeNormals,
                const labelList& meshEdges,
                labelList& pointType,
                labelList& edgeType,
                labelList& faceType,
                DynamicList<label>& boundaryCornerPoints
            );

            label modifyBoundaryFaces
            (
                const indirectPrimitivePatch& pp,
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& faceType,
                const List<labelPair>& newFaceCells,
                const List<labelPair>& newEdgeCells,
                const List<labelPair>& newPointCells,
                polyTopoChange& meshMod
            );

            label addTopFaces
            (
                const indirectPrimitivePatch& pp,
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& faceType,
                const List<labelPair>& newFaceCells,
                const List<labelPair>& newEdgeCells,
                const List<labelPair>& newPointCells,
                const labelList& addedPoints,
                const List<DynamicList<labelPair>>& convexEdgePts,
                const List<DynamicList<labelPair>>& convexFacePts,
                const labelList& baffleLevel,
                polyTopoChange& meshMod
            );

            label addSideFaces
            (
                const indirectPrimitivePatch& pp,
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& faceType,
                const List<labelPair>& newFaceCells,
                const List<labelPair>& newEdgeCells,
                const List<labelPair>& newPointCells,
                const labelList& addedPoints,
                const List<DynamicList<labelPair>>& convexEdgePts,
                const List<DynamicList<labelPair>>& convexFacePts,
                const labelList& baffleLevel,
                const labelList& meshEdges,
                polyTopoChange& meshMod
            );

            label modifyOtherFaces
            (
                const labelList& pointType,
                const labelList& edgeType,
                const labelList& faceType,
                const labelList& addedPoints,
                const List<DynamicList<labelPair>>& convexEdgePts,
                const List<DynamicList<labelPair>>& convexFacePts,
                const labelList& baffleLevel,
                polyTopoChange& meshMod
            );

            List<labelPair> createBaffles();

            void mergeRefinementBaffles();

            void setRefinement();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

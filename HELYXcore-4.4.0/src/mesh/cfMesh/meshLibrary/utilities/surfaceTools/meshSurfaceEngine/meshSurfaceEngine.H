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
    meshSurfaceEngine

Description
    Calculates surface of the mesh

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngine_H
#define meshSurfaceEngine_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/Time/Time.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "containers/Lists/SubList/SubList.H"
#include "primitives/bools/lists/boolList.H"
#include "utilities/containers/DynList/DynList.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "meshes/primitiveShapes/triangle/triangle.H"
#include "containers/HashTables/Map/Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meshSurfaceEngine Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngine
{
    // Private data

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- Number of active patch
        const label activePatch_;

        //- Boundary points
        mutable labelList* boundaryPointsPtr_;

        //- Boundary faces
        mutable faceList::subList* boundaryFacesPtr_;

        //- Patches boundary faces are in
        mutable labelList* boundaryFacePatchPtr_;

        //- Face owners
        mutable labelList* boundaryFaceOwnersPtr_;

        //- Point faces addressing
        mutable VRWGraph* pointFacesPtr_;
        mutable VRWGraph* pointInFacePtr_;

        //- Point-patches addressing
        mutable VRWGraph* pointPatchesPtr_;

        //- PointBoundaryPoint addressing
        mutable labelList* bppPtr_;

        //- Point points addressing
        mutable VRWGraph* pointPointsPtr_;

        //- Edges
        mutable edgeList* edgesPtr_;

        //- Boundary point-edges addressing
        mutable VRWGraph* bpEdgesPtr_;

        //- Edge faces addressing
        mutable VRWGraph* edgeFacesPtr_;

        //- Face edges addressing
        mutable VRWGraph* faceEdgesPtr_;

        //- Edge-patches addressing
        mutable VRWGraph* edgePatchesPtr_;

        //- Face-faces addressing
        mutable VRWGraph* faceFacesPtr_;

        //- Point normals
        mutable vectorField* pointNormalsPtr_;

        //- Face normals
        mutable vectorField* faceNormalsPtr_;

        //- Face centres
        mutable vectorField* faceCentresPtr_;


    // Private data for parallel execution

        //- Global boundary point label
        mutable labelList* globalBoundaryPointLabelPtr_;

        //- Global boundary point to local addressing
        mutable Map<label>* globalBoundaryPointToLocalPtr_;

        //- Boundary point-processors addressing
        mutable VRWGraph* bpProcsPtr_;

        //- Neighbour processors for communication when sending point data
        mutable DynList<label>* bpNeiProcsPtr_;

        //- Global boundary edge label
        mutable labelList* globalBoundaryEdgeLabelPtr_;

        //- Global boundary edge to local addressing
        mutable Map<label>* globalBoundaryEdgeToLocalPtr_;

        //- Boundary edge-processors addressing
        mutable VRWGraph* beProcsPtr_;

        //- Neighbour processors for communication when sending edge data
        mutable DynList<label>* beNeiProcsPtr_;

        //- Processor containing other face and face-patch addressing
        mutable Map<label>* otherEdgeFaceAtProcPtr_;
        mutable Map<label>* otherEdgeFacePatchPtr_;

        //- Global label for boundary faces
        mutable labelList* globalBoundaryFaceLabelPtr_;


    // Private member functions

        //- calculate boundary nodes, faces and addressing
        void calculateBoundaryFaces() const;
        void calculateBoundaryOwners() const;
        void calculateBoundaryNodes() const;

        void calculateBoundaryFacePatches() const;

        void calculatePointFaces() const;

        void calculatePointPatches() const;

        void calculatePointPoints() const;

        void calculatePointNormals() const;
        void calculateFaceNormals() const;
        void calculateFaceCentres() const;

        void updatePointNormalsAtProcBoundaries() const;

        void calculateEdgesAndAddressing() const;
        void calculateFaceEdgesAddressing() const;
        void calculateEdgeFacesAddressing() const;
        void calculateEdgePatchesAddressing() const;
        void calculateFaceFacesAddressing() const;


    // Private member functions for parallel runs

        void calcGlobalBoundaryPointLabels() const;
        void calcGlobalBoundaryEdgeLabels() const;
        void calcAddressingForProcEdges() const;
        void calcGlobalBoundaryFaceLabels() const;


    // Private member functions

        void clearOut();


public:

        // Declare name of the class and its debug switch
        ClassName("meshObject");

        friend class meshSurfaceEngineModifier;


    // Constructors

        //- Construct from mesh
        meshSurfaceEngine(polyMeshGen& mesh);

        //- Construct from mesh and patch
        meshSurfaceEngine(polyMeshGen& mesh, const label patchI);

        //- Disallow default bitwise copy construction
        meshSurfaceEngine(const meshSurfaceEngine&) = delete;


    //- Destructor
    ~meshSurfaceEngine();


    // Member Functions

        inline const polyMeshGen& mesh() const;
        inline const pointFieldPMG& points() const;
        inline const faceListPMG& faces() const;
        inline const cellListPMG& cells() const;

        inline const labelList& boundaryPoints() const;

        inline const faceList::subList& boundaryFaces() const;

        //- Patch label for each boundary face
        inline const labelList& boundaryFacePatches() const;

        inline const labelList& faceOwners() const;

        inline const VRWGraph& pointFaces() const;
        inline const VRWGraph& pointInFaces() const;

        inline const VRWGraph& pointPoints() const;

        inline const vectorField& pointNormals() const;

        inline const vectorField& faceNormals() const;

        inline const vectorField& faceCentres() const;

        inline const labelList& bp() const;

        inline const edgeList& edges() const;

        inline const VRWGraph& boundaryPointEdges() const;

        inline const VRWGraph& edgeFaces() const;

        inline const VRWGraph& faceEdges() const;

        inline const VRWGraph& edgePatches() const;

        inline const VRWGraph& faceFaces() const;


    // Functions for parallel runs

        //- Global boundary point label
        inline const labelList& globalBoundaryPointLabel() const;

        //- Global point label to local label. Only for processors points
        inline const Map<label>& globalToLocalBndPointAddressing() const;

        //- Processors which contain the vertex
        inline const VRWGraph& bpAtProcs() const;

        //- Communication matrix for sending point data
        inline const DynList<label>& bpNeiProcs() const;

        //- Global boundary edge label
        inline const labelList& globalBoundaryEdgeLabel() const;

        //- Global boundary edge label to local label. Only for processor edges
        inline const Map<label>& globalToLocalBndEdgeAddressing() const;

        //- Processors which contain the edges
        inline const VRWGraph& beAtProcs() const;

        //- Communication matrix for sending edge data
        inline const DynList<label>& beNeiProcs() const;

        //- Patch label and processor label at which the other surface face
        //  sharing processor edge is located. Only for processor edges
        inline const Map<label>& otherEdgeFaceAtProc() const;
        inline const Map<label>& otherEdgeFacePatch() const;

        //- Global boundary face label
        inline const labelList& globalBoundaryFaceLabel() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngine&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

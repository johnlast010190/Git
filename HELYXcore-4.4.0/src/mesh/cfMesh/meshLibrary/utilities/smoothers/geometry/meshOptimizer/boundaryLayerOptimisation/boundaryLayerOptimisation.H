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
    bondaryLayerOptimisation

Description
    A class intended for improving quality of boundary layer cells. It is intended
    to be used on meshes with one boundary layer which can later be refined.

SourceFiles
    boundaryLayerOptimisation.C
    boundaryLayerOptimisationFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryLayerOptimisation_H
#define boundaryLayerOptimisation_H

#include "utilities/containers/DynList/DynList.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"
#include "utilities/containers/Lists/edgeLongList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class dictionary;

/*---------------------------------------------------------------------------*\
                Class boundaryLayerOptimisation Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayerOptimisation
{
    // Private data

        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- Const pointer to meshSurfaceEngine
        mutable const meshSurfaceEngine* meshSurfacePtr_;
        const bool deleteMeshSurface_;

        //- Mesh surface partitioner
        mutable const meshSurfacePartitioner* partitionerPtr_;

        //- Boundary layer hairs
        edgeLongList hairEdges_;

        //- Hair edges attached to a boundary point
        VRWGraph hairEdgesAtBndPoint_;

        //- Hair edge to other hair edges
        VRWGraph hairEdgesNearHairEdge_;

        //- Is boundary face a base for a prism in the bnd layer
        boolList isBndLayerBase_;

        //- Is boundary face part of a layer where a layer exits
        boolList isExitFace_;

        //- Classification of hair edges
        List<direction> hairEdgeType_;

        //- Stores information where boundary hairs are made thinner
        boolList thinnedHairEdge_;

        //- Maximum number of iterations
        label maxNumIterations_;

        //- Number of iterations for smoothing of hairs
        label nSmoothNormals_;

        //- Relative tolerance for thickness variation
        scalar relThicknessTol_;

        //- Feature size factor used for curvature-based optimisation
        scalar featureSizeFactor_;

        //- Activate calculation of normals
        bool reCalculateNormals_;


    // Private member functions

        //- Access to mesh surface
        const meshSurfaceEngine& meshSurface() const;

        //- Access to meshSurfacePartitioner
        const meshSurfacePartitioner& surfacePartitioner() const;

        //- Calculate normal vector for all patches at a point
        //  points can be filtered by a type of hair edges
        typedef std::map<label, std::pair<point, scalar>> patchNormalType;
        typedef std::map<label, patchNormalType> pointNormalsType;
        void calculateNormalVectors
        (
            const direction eType,
            pointNormalsType&
        ) const;

        //- Calculate normal vectors
        void calculateNormalVectorsSmother
        (
            const direction eType,
            pointNormalsType&
        );

        //- Calculate hairEdges
        void calculateHairEdges();

        //- Calculate hair vectors at the boundary
        void calculateHairVectorsAtTheBoundary(vectorField&);

        //- Optimise hair normals at the boundary
        //  these are the regions where boundary layers exit the domain
        void optimiseHairNormalsAtTheBoundary();

        //- Optimise hair normals inside the mesh
        void optimiseHairNormalsInside();

        //- Calculate hair edges at a boundary faces
        void hairEdgesAtBndFace
        (
            const label cellI,
            const label baseFaceI,
            DynList<edge>&
        ) const;

        //- Calculate the thickness of the layer for edge heI
        //  respective to the height of edge heJ
        scalar calculateThickness
        (
            const label heI,
            const label heJ
        ) const;

        //- Calculate thickness of the layer by checking intersections
        //  between the hair edge and the quad faces extruded from
        //  nighbouring faces
        scalar calculateThicknessOverCell
        (
            const label heI,
            const label cellI,
            const label baseFaceI
        ) const;

        //- Optimise thickness variation
        void optimiseThicknessVariation
        (
            const direction edgeType = (INSIDE|BOUNDARY)
        );

        //- Optimise layers at exitting faces due to requests from the inside
        bool optimiseLayersAtExittingFaces();

        //- Write vectors into a VTK file. Helper for debugging
        static void writeVTK
        (
            const fileName& fName,
            const pointField& origin,
            const vectorField& vecs
        );

        //- Write vector correcposing to hair edges. Helper for debugging
        void writeHairEdges
        (
            const fileName& fName,
            const direction eType,
            const vectorField& vecs
        ) const;

        //- Write selected hair edges into a file
        void writeHairEdges(const fileName& fName, const direction eType) const;


    // Private enumerators

        enum hairEdgeTypes_
        {
            NONE = 0,
            ATEDGE = 1,
            ATCORNER = 2,
            BOUNDARY = 4,
            INSIDE = 8,
            FEATUREEDGE = 16
        };


public:

    // Constructors

        //- Construct from polyMeshGen
        boundaryLayerOptimisation(polyMeshGen& mesh);

        //- Construct from polyMeshGen and meshSurfaceEngine
        boundaryLayerOptimisation
        (
            polyMeshGen& mesh,
            const meshSurfaceEngine&
        );

        //- Disallow default bitwise copy construction
        boundaryLayerOptimisation(const boundaryLayerOptimisation&) = delete;


    //- Destructor
    ~boundaryLayerOptimisation();


    // Member Functions

        //- Set the maximum number of iterations
        void setMaxNumIterations(const label maxNumIterations);

        //- Set the number of normal smoothing iterations (default is 5)
        void setNumNormalsSmoothingIterations(const label nSmoothNormal);

        //- Shall normals be re-calculated (default true)
        void recalculateNormals(const bool);

        //- Set the relative thickness tolerance (default 0.15)
        void setRelativeThicknessTolerance(const scalar);

        //- Set the feature size factor (default 0.3)
        void setFeatureSizeFactor(const scalar);

        //- Return hair edges
        const edgeLongList& hairEdges() const;

        //- Hair edges attached to a boundary point
        const VRWGraph& hairEdgesAtBndPoint() const;

        //- Boundary faces which serve as base faces of a boundary layer
        //  are set to true
        const boolList& isBaseFace() const;

        //- Boundary faces where the layers exit at the boundary
        const boolList& isExitFace() const;

        //- Performs boundary layer optimisation
        void optimiseLayer();


    // Static member functions

        //- Read the settings from dictionary
        static void readSettings
        (
            const dictionary&,
            boundaryLayerOptimisation&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const boundaryLayerOptimisation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

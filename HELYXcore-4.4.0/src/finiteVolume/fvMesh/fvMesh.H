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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2011-2024 OpenFOAM Foundation
    (c) 2010-2025 Engys Ltd.

Class
    Foam::fvMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

    NOTE ON USAGE:
    fvMesh contains all the topological and geometric information
    related to the mesh.  It is also responsible for keeping the data
    up-to-date.  This is done by deleting the cell volume, face area,
    cell/face centre, addressing and other derived information as
    required and recalculating it as necessary.  The fvMesh therefore
    reserves the right to delete the derived information upon every
    topological (mesh refinement/morphing) or geometric change (mesh
    motion).  It is therefore unsafe to keep local references to the
    derived data outside of the time loop.

SourceFiles
    fvMesh.C
    fvMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fvMesh_H
#define fvMesh_H

#include "meshes/polyMesh/polyMesh.H"
#include "meshes/lduMesh/lduMesh.H"
#include "meshes/data/data.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolation.H"
#include "finiteVolume/fvSchemes/fvSchemes.H"
#include "finiteVolume/fvSolution/fvSolution.H"
#include "fields/DimensionedFields/SlicedDimensionedField/SlicedDimensionedField.H"
#include "fields/volFields/slicedVolFieldsFwd.H"
#include "fields/surfaceFields/slicedSurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshLduAddressing;
class fvMeshStitcher;
class fvMeshTopoChanger;
class fvMeshDistributor;
class fvMeshMover;
class fvMeshGIBChanger;
class volMesh;
class polyDistributionMap;
template<class Type>
class fvMatrix;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricBoundaryField;

template<class T, class Container>
class CompactListList;

namespace fvMeshTopoChangers
{
    class list;
}

/*---------------------------------------------------------------------------*\
                            Class fvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMesh
:
    public polyMesh,
    public lduMesh,
    public surfaceInterpolation,
    public data
{
public:

    // Public data types

        //- Extent to which to stitch on read and readUpdate. By default, a
        //  full geometric stitch is performed. A non-geometric stitch can be
        //  done as an optimisation in situations when finite volume geometry
        //  is not needed (e.g. decomposition). Stitching can also be
        //  prevented altogether if that is appropriate (e.g. if the mesh is
        //  loaded for mapping in an un-stitched state).
        enum class stitchType
        {
            none,
            nonGeometric,
            geometric
        };


private:

    // Private data

        //- Boundary mesh
        fvBoundaryMesh boundary_;

        //- The stitcher function class
        autoPtr<fvMeshStitcher> stitcher_;

        //- The topo-changer function class
        autoPtr<fvMeshTopoChanger> topoChanger_;

        //- The distributor function class
        autoPtr<fvMeshDistributor> distributor_;

        //- The mover function class
        autoPtr<fvMeshMover> mover_;

        //- The GIB-changer function class
        autoPtr<fvMeshGIBChanger> GIBChanger_;


    // Demand-driven data

        //- Matrix addressing
        mutable fvMeshLduAddressing* lduPtr_;

        //- Face-poly-face addressing
        mutable GeometricBoundaryField<label, fvsPatchField, surfaceMesh>*
            polyFacesBfPtr_;

        //- Offsets for poly-bFace-patch and patch-face addressing
        mutable labelList* polyBFaceOffsetsPtr_;

        //- Poly-bFace-patch addressing
        mutable labelList* polyBFaceOffsetPatchesPtr_;

        //- Poly-bFace-patch-face addressing
        mutable labelList* polyBFaceOffsetPatchFacesPtr_;

        //- Poly-bFace-patch addressing
        mutable CompactListList<label, labelList>* polyBFacePatchesPtr_;

        //- Poly-bFace-patch-face addressing
        mutable CompactListList<label, labelList>* polyBFacePatchFacesPtr_;

        //- Face-owner addressing
        mutable GeometricBoundaryField<label, fvsPatchField, surfaceMesh>*
            ownerBfPtr_;

        //- Current time index for cell volumes
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Cell volumes
        mutable SlicedDimensionedField<scalar, volMesh>* VPtr_;

        //- Cell volumes old time level
        mutable DimensionedField<scalar, volMesh>* V0Ptr_;

        //- Cell volumes old-old time level
        mutable DimensionedField<scalar, volMesh>* V00Ptr_;

        //- Face area vectors
        mutable slicedSurfaceVectorField* SfSlicePtr_;

        //- Face area vectors
        mutable surfaceVectorField* SfPtr_;

        //- Mag face area vectors
        mutable slicedSurfaceScalarField* magSfSlicePtr_;

        //- Mag face area vectors
        mutable surfaceScalarField* magSfPtr_;

        //- Cell centres
        mutable slicedVolVectorField* CSlicePtr_;

        //- Cell centres
        mutable volVectorField* CPtr_;

        //- Face centres
        mutable slicedSurfaceVectorField* CfSlicePtr_;

        //- Face centres
        mutable surfaceVectorField* CfPtr_;

        //- Face motion fluxes
        mutable surfaceScalarField* phiPtr_;

        //- fvSchemes created on demand
        mutable autoPtr<fvSchemes> fvSchemes_;

        //- fvSolution created on demand
        mutable autoPtr<fvSolution> fvSolution_;



    // Private Member Functions

        // Storage management

            //- Clear geometry but not the old-time cell volumes
            void clearGeomNotOldVol();

            //- Clear geometry like clearGeomNotOldVol but recreate any
            //  geometric demand-driven data that was set
            void updateGeomNotOldVol();

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Store old times fields for all fields of a given type and class
            template<class Type, template<class> class GeoField>
            void storeOldTimeFields();

            //- Store old times fields for all fields of a given class
            template<template<class> class GeoField>
            void storeOldTimeFields();

            //- Store old times fields for all fields.
            //  Used prior to field mapping following mesh-motion
            //  and NCC stitching.
            void storeOldTimeFields();

            //- Null oldest time field for all fields of a given type and class
            template<class Type, template<class> class GeoField>
            void nullOldestTimeFields();

            //- Null oldest time field for all fields of a given class
            template<template<class> class GeoField>
            void nullOldestTimeFields();

            //- Null oldest time field for all fields.
            //  Used prior to field mapping following mesh topology change
            void nullOldestTimeFields();


        // Make topological data

            //- Get the patch types for the poly faces field
            wordList polyFacesPatchTypes() const;

            //- Modify face-poly-face addressing
            GeometricBoundaryField<label, fvsPatchField, surfaceMesh>&
                polyFacesBfRef();


        // Make geometric data

            //- Make the face area vectors
            void makeSf() const;

            //- Make the mag face area vectors
            void makeMagSf() const;

            //- Make the cell centres
            void makeC() const;

            //- Make the face centres
            void makeCf() const;

            //- Modify cell face area vectors
            surfaceVectorField& SfRef();

            //- Modify cell face area magnitudes
            surfaceScalarField& magSfRef();

            //- Modify cell centres
            volVectorField& CRef();

            //- Modify face centres
            surfaceVectorField& CfRef();


public:

    // Public typedefs

        typedef fvMesh Mesh;
        typedef fvBoundaryMesh BoundaryMesh;

        friend fvMeshTopoChangers::list;


    // Declare name of the class and its debug switch
    ClassName("fvMesh");


    // Static data

        //- Set of names of registered geometric fields:
        //  V, Sf, magSf, C, Cf
        const static HashSet<word> geometryFields;


    // Constructors

        //- Construct from IOobject
        //  with the option to instantiate the mesh changers and the option
        //  to prevent some or all of the stitching
        explicit fvMesh
        (
            const IOobject& io,
            const bool changers = false,
            const stitchType stitch = stitchType::geometric
        );

        //- Move constructor from components without boundary.
        //  Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Move constructor without boundary from cells rather than
        //  owner/neighbour. Boundary is added using addPatches()
        //  member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true,
            const bool autoWrite = true
        );

        //- Move constructor from blockMesh
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const PtrList<dictionary>& boundaryDicts,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const bool syncPar = true
        );

        //- Move constructor
        fvMesh(fvMesh&&);

        //- Disallow default bitwise copy construction
        fvMesh(const fvMesh&) = delete;


    //- Destructor
    virtual ~fvMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addFvPatches
            (
                const List<polyPatch*>&,
                const bool validBoundary = true
            );

            //- Update the mesh based on the mesh files saved in time
            //  directories
            readUpdateState readUpdate
            (
                const stitchType stitch = stitchType::geometric
            );


        // Access

            //- Return the top-level database
            const Time& time() const
            {
                return polyMesh::time();
            }

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const
            {
                return true;
            }

            //- Return the object registry - resolve conflict polyMesh/lduMesh
            virtual const objectRegistry& thisDb() const
            {
                return polyMesh::thisDb();
            }

            //- Return reference to name
            //  Note: name() is currently ambiguous due to derivation from
            //  surfaceInterpolation
            const word& name() const
            {
                return polyMesh::name();
            }

            //- Return reference to boundary mesh
            const fvBoundaryMesh& boundary() const
            {
                return boundary_;
            }

            //- Return reference to polyMesh
            const polyMesh& operator()() const
            {
                return *this;
            }

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return polyMesh::comm();
            }


            // Overlap

                //- Interpolate interpolationCells only
                virtual void interpolate(volScalarField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volVectorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volSphericalTensorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volSymmTensorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(scalarField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(vectorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(sphericalTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(symmTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(tensorField&) const
                {}

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<scalar> solve
                (
                    fvMatrix<scalar>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<vector> solve
                (
                    fvMatrix<vector>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<sphericalTensor> solve
                (
                    fvMatrix<sphericalTensor>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<symmTensor> solve
                (
                    fvMatrix<symmTensor>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<tensor> solve
                (
                    fvMatrix<tensor>&,
                    const dictionary&
                ) const;


            //- Internal face owner. Note bypassing virtual mechanism so
            //  e.g. relaxation always gets done using original addressing
            const labelUList& owner() const
            {
                return fvMesh::lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return fvMesh::lduAddr().upperAddr();
            }

            //- Return whether the fvMesh is conformal with the polyMesh
            bool conformal() const
            {
                return !(polyFacesBfPtr_ && SfPtr_);
            }

            //- Get the poly faces IO object
            IOobject polyFacesBfIO(IOobject::readOption r) const;

            //- Return face-poly-face addressing
            const GeometricBoundaryField<label, fvsPatchField, surfaceMesh>&
                polyFacesBf() const;

            //- Return poly-bFace-patch addressing
            const CompactListList<label, labelList>& polyBFacePatches() const;

            //- Return poly-bFace-patch-face addressing
            const CompactListList<label, labelList>&
                polyBFacePatchFaces() const;

            //- Return face-owner addressing
            const GeometricBoundaryField<label, fvsPatchField, surfaceMesh>&
                ownerBf() const;

            //- Return the stitcher function class
            const fvMeshStitcher& stitcher() const;

            //- Return the topo-changer function class
            const fvMeshTopoChanger& topoChanger() const;

            //- Return the distributor function class
            const fvMeshDistributor& distributor() const;

            //- Return the mover function class
            const fvMeshMover& mover() const;

            //- Return the GIB-changer function class
            const fvMeshGIBChanger& GIBChanger() const;

            //- Return cell volumes
            virtual const DimensionedField<scalar, volMesh>& V() const;

            //- Return old-time cell volumes
            virtual const DimensionedField<scalar, volMesh>& V0() const;

            //- Return old-old-time cell volumes
            virtual const DimensionedField<scalar, volMesh>& V00() const;

            //- Return sub-cycle cell volumes
            virtual tmp<DimensionedField<scalar, volMesh>> Vsc() const;

            //- Return sub-cycl old-time cell volumes
            virtual tmp<DimensionedField<scalar, volMesh>> Vsc0() const;

            //- Return cell face area vectors
            virtual const surfaceVectorField& Sf() const;

            //- Return cell face area magnitudes
            virtual const surfaceScalarField& magSf() const;

            //- Return cell face motion fluxes
            const surfaceScalarField& phi() const;

            //- Return cell centres
            virtual const volVectorField& C() const;

            //- Return face centres
            virtual const surfaceVectorField& Cf() const;

            //- Return the list of fields of type GeoField,
            //  i.e. fields that are not geometry.
            template<class GeoField>
            UPtrList<GeoField> fields(const bool strict = false) const;

            //- Return the list of current fields of type GeoField,
            //  i.e. fields that are not old-times or geometry.
            template<class GeoField>
            UPtrList<GeoField> curFields() const;

            //- Return face deltas as surfaceVectorField
            tmp<surfaceVectorField> delta() const;

            //- Return a labelType of valid component indicators
            //  1 : valid (solved)
            // -1 : invalid (not solved)
            // e.g. for tensor with z not solved: ( 1 1 -1 1 1 -1 -1 -1 1)
            template<class Type>
            typename pTraits<Type>::labelType validComponents() const;

            //- Return a labelType of valid component indicators
            //  1 : valid (solved)
            //  0 : invalid (not solved)
            // e.g. for tensor with z not solved: ( 1 1 0 1 1 0 0 0 0)
            template<class Type>
            typename pTraits<Type>::labelType validComponents2() const;

            //- Make tensor an identity in the empty directions so that inverse
            //  can be taken
            template<class Type>
            void stabiliseEmptyDirections(Field<Type>& t) const;

            template<class Type, template<class> class PatchField, class GeoMesh>
            void stabiliseEmptyDirections
            (
                GeometricField<Type, PatchField, GeoMesh>& t
            ) const;


            //- Find patch indices for a given active polyPatch type
            template<class Type>
            labelHashSet findActivePatchIDs() const;

            //- Return the fvSchemes
            const fvSchemes& schemes() const;

            //- Return the fvSchemes
            const fvSolution& solution() const;

            //- Return the fvSchemes
            fvSchemes& schemes();

            //- Return the fvSchemes
            fvSolution& solution();


        // Edit

            //- Is the mesh instantiated with the new 'changers' structure
            bool hasChangers() const;

            //- Is the construction of non-conformal couples needed
            bool needsStitching() const;

            //- Does the mesh topology change?
            bool topoChanging() const;

            //- Is this mesh dynamic?
            bool dynamic() const;

            //- Update the mesh for both topology change and mesh motion
            virtual bool update();

            //- Update the mesh for topology change only
            virtual bool topoChange();

            //- Move the mesh
            virtual bool move();

            //- Clear mesh motion flux
            void clearMeshPhi();
            void clearV0();
            void clearGeomData();

            //- Clear all geometry and addressing
            void clearOut();

            //- Update mesh corresponding to the given map
            virtual void topoChange(const polyTopoChangeMap& map);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& mdpm);

            //- Redistribute or update using the given map
            virtual void updateGIB();

            //- Redistribute or update using the given map
            virtual void storeGIBFields();

            //- Reset the points without storing old points
            //  or returning swept volumes
            virtual void setPoints(const pointField&);

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const pointField&);

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> moveGIBPoints(const pointField&);

            //- Return old points for mesh motion
            virtual const pointField& oldPoints() const;

            //- Conform the fvMesh to the polyMesh
            void conform
            (
                const surfaceScalarField& phi =
                    NullObjectRef<surfaceScalarField>()
            );

            //- Unconform the fvMesh from the polyMesh
            void unconform
            (
                const GeometricBoundaryField
                <
                    label,
                    fvsPatchField,
                    surfaceMesh
                >& polyFacesBf,
                const surfaceVectorField& Sf,
                const surfaceVectorField& Cf,
                const surfaceScalarField& phi =
                    NullObjectRef<surfaceScalarField>(),
                const bool sync = true
            );

            //- Map all fields in time using given map.
            void mapFields(const polyTopoChangeMap& map);

            //- Add/insert single patch. If validBoundary the new situation
            //  is consistent across processors.
            virtual void addPatch
            (
                const label insertPatchi,
                const polyPatch& patch,
                const dictionary& patchFieldDict,
                const word& defaultPatchFieldType,
                const bool validBoundary
            );

            //- Reorder and trim existing patches. If validBoundary the new
            //  situation is consistent across processors
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Remove boundary patches. Warning: fvPatchFields hold ref to
            //  these fvPatches.
            void removeFvBoundary();

            //- Swap mesh
            //  For run-time mesh replacement and mesh to mesh mapping
            void swap(fvMesh&);

            //- Modify cell face motion fluxes
            surfaceScalarField& setPhi();

            //- Modify old-time cell volumes
            DimensionedField<scalar, volMesh>& setV0();

            //- Correct velocity for GIB
            virtual tmp<vectorField> velocityCorrect(const vectorField&) const;


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write
            ) const;

            //- Write mesh using IO settings from time
            virtual bool write(const bool write = true) const;


    // Member Operators

        bool operator!=(const fvMesh&) const;
        bool operator==(const fvMesh&) const;

        //- Disallow default bitwise assignment
        void operator=(const fvMesh&) = delete;
};


template<>
//typename pTraits<sphericalTensor>::labelType
pTraits<sphericalTensor>::labelType
fvMesh::validComponents<sphericalTensor>() const;


template<>
pTraits<sphericalTensor>::labelType
fvMesh::validComponents2<sphericalTensor>() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/fvMeshTemplates.C"
    #include "fvMesh/fvPatches/fvPatch/fvPatchFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

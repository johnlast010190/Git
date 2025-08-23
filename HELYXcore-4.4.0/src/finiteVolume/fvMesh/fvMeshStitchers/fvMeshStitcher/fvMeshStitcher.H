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
    (c) 2021-2023 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::fvMeshStitcher

Description
    Mesh manipulator that uses the intersection provided by the cyclic
    non-conformal poly patches to create non-conformal finite volume
    interfaces.

SourceFiles
    fvMeshStitcher.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshStitcher_H
#define fvMeshStitcher_H

#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/Fields/fieldMappers/setSizeFieldMapper/setSizeFieldMapper.H"
#include "patchToPatch/intersection/intersectionPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCyclicFvPatch;

/*---------------------------------------------------------------------------*\
                        Class fvMeshStitcher Declaration
\*---------------------------------------------------------------------------*/

class fvMeshStitcher
{
protected:

    // Protected Typedefs

        //- Alias for surface boundary fields
        template<class Type>
        using SurfaceFieldBoundary =
            GeometricBoundaryField<Type, fvsPatchField, surfaceMesh>;

        //- Alias the patchToPatch intersection part struct
        using part = patchToPatches::intersection::part;

        //- Alias the patchToPatch intersection couple struct
        using couple = patchToPatches::intersection::couple;


    // Protected Classes

        //- Set size fvPatchFieldMapper
        struct fvPatchFieldSetSizer
        :
            public fvPatchFieldMapper,
            public setSizeFieldMapper
        {
            fvPatchFieldSetSizer(const label size)
            :
                setSizeFieldMapper(size)
            {}
        };


private:

    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;


    // Private Member Functions

        // Intersection

            //- Perform intersections for the given non-conformal interface
            void intersectNonConformalCyclic
            (
                const nonConformalCyclicFvPatch& nccFvp,
                SurfaceFieldBoundary<label>& polyFacesBf,
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                const tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf,
                List<part>& origEdgeParts
            ) const;

            //- Convert per-patch lists of edge parts into a single
            //  owner-orig-boundary list of edge parts that can be applied to
            //  the mesh
            List<part> calculateOwnerOrigBoundaryEdgeParts
            (
                const List<List<part>>& patchEdgeParts
            ) const;

            //- Apply a list of owner-orig-boundary edge parts to
            //  non-orig-patch faces (both internal and boundary)
            void applyOwnerOrigBoundaryEdgeParts
            (
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                const List<part>& ownerOrigBoundaryEdgeParts
            ) const;

            //- Stabilise the faces that have had non-conformal coupled parts
            //  cut out from them
            void stabiliseOrigPatchFaces
            (
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf
            ) const;

            //- Perform intersections for all non-conformal interfaces
            void intersectNonConformalCyclics
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                const bool haveTopology
            ) const;

            //- Apply the given non-conformal topology (polyFacesBf) to the
            //  face geometry fields for the given patch type. Just set small
            //  stabilisation values. Don't actually do any cutting. Used for
            //  initialisation, or in situations when cutting cannot be done
            //  (e.g., when post processing a processor case).
            template<class NonConformalFvPatch>
            inline void createNonConformalStabilisationGeometry
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            ) const;

            //- Perform intersections for non-conformal discrete mixing
            //  interfaces
            void intersectNonConformalDiscreteMixing
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            );


        // Motion

            //- Return whether or not this stitcher supports a changing mesh
            virtual bool changing() const = 0;

            //- Correct the mesh fluxes following a conform operation
            virtual void conformCorrectMeshPhi(surfaceScalarField& phi) = 0;

            //- Initialise correction of the mesh fluxes
            virtual void createNonConformalCorrectMeshPhiGeometry
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            ) = 0;

            //- Correct the mesh fluxes following an unconform operation
            virtual void unconformCorrectMeshPhi
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                surfaceScalarField& phi
            ) = 0;


        // Field Mapping

            //- Resize the patch fields of a given type and class to match the
            //  sizes of the patches in the mesh
            template<class Type, template<class> class GeoField>
            void resizePatchFields();

            //- As above, for all types.
            template<template<class> class GeoField>
            void resizePatchFields();

            //- Pre-conform surface fields of a given type by separating NCC
            //  and original parts of non-conformal couplings and storing them
            //  in the database as separate fields
            template<class Type>
            void preConformSurfaceFields();

            //- As above, for all types.
            void preConformSurfaceFields();

            //- Post-unconform surface fields of a given type by looking up
            //  NCC and original parts of non-conformal couplings and combining
            //  them into a single non-conformal boundary field
            template<class Type>
            void postUnconformSurfaceFields();

            //- As above, for all types.
            void postUnconformSurfaceFields();

            //- Evaluate all non-conformal vol patch fields of a given type
            template<class Type>
            void evaluateVolFields();

            //- As above, for all types.
            void evaluateVolFields();

            //- Special post-unconform for surface velocities
            void postUnconformSurfaceVelocities();


protected:

    // Protected Member Functions

        //- Access the boundary field reference of a field, without updating
        //  the time index, storing old time fields etc.
        template<class GeoField>
        static typename GeoField::Boundary& boundaryFieldRefNoUpdate
        (
            GeoField& fld
        );


        // Checking

            //- Return the non-dimensional cell openness for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> openness() const;

            //- Return the non-dimensional old-time volume conservation error
            //  for a specified old-time index for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> volumeConservationError
            (
                const label n
            ) const;


public:

    //- Runtime type information
    TypeName("fvMeshStitcher");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshStitcher,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshStitcher(fvMesh&);

        //- Disallow default bitwise copy construction
        fvMeshStitcher(const fvMeshStitcher&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshStitcher
        static autoPtr<fvMeshStitcher> New(fvMesh&, const bool changers);


    //- Destructor
    virtual ~fvMeshStitcher();


    // Member Functions

        //- Return the fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Does this stitcher do anything?
        bool stitches() const;

        //- Disconnect the mesh by removing faces from the nonConformalCyclics
        bool disconnect(const bool changing, const bool geometric);

        //- Connect the mesh by adding faces into the nonConformalCyclics
        bool connect
        (
            const bool changing,
            const bool geometric,
            const bool load
        );

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}

        //- Write the stitcher state
        virtual bool write(const bool write = true) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshStitcher&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcherMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

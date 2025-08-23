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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::meshAdder

Description
    Adds two fvMeshes without using any polyMesh morphing.
    Uses polyMeshAdder.

SourceFiles
    meshAdder.C
    meshAdderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef meshAdder_H
#define meshAdder_H

#include "polyMeshAdder/polyMeshAdder.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldsFwd.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFieldsFwd.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class IOobject;
class faceCoupleInfo;
class IOobjectList;
class fvMesh;
class volMesh;
class surfaceMesh;
class mapAddedPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class meshAdder Declaration
\*---------------------------------------------------------------------------*/

class meshAdder
:
    public polyMeshAdder
{

private:

    // Private Member Functions

        //- Calculate map from new patch faces to old patch faces. -1 where
        //  could not map.
        static labelList calcPatchMap
        (
            const label oldStart,
            const label oldSize,
            const labelList& oldToNew,
            const polyPatch& newPatch,
            const label unmappedIndex
        );

        //- Update single volField.
        template<class Type>
        static void MapVolField
        (
            const mapAddedPolyMesh& meshMap,

            VolField<Type>& fld,
            const VolField<Type>& fldToAdd,
            const bool fullyMapped
        );

        //- Update single surfaceField.
        template<class Type>
        static void MapSurfaceField
        (
            const mapAddedPolyMesh& meshMap,
            SurfaceField<Type>& fld,
            const SurfaceField<Type>& fldToAdd,
            const bool fullyMapped
        );

        //- Update single pointField.
        template<class Type>
        static void MapPointField
        (
            const pointMesh& mesh,
            const mapAddedPolyMesh& meshMap,
            const labelListList& oldMeshPoints,

            GeometricField<Type, pointPatchField, pointMesh>& fld,
            const GeometricField<Type, pointPatchField, pointMesh>& fldToAdd
        );

        //- Update single dimensionedField.
        template<class Type>
        static void MapDimField
        (
            const mapAddedPolyMesh& meshMap,
            DimensionedField<Type, volMesh>& fld,
            const DimensionedField<Type, volMesh>& fldToAdd
        );

public:

    // Declare name of the class and its debug switch
    ClassName("meshAdder");


    // Member Functions

        //- Inplace add mesh to fvMesh. Maps all stored fields. Returns map.
        static autoPtr<mapAddedPolyMesh> add
        (
            fvMesh& mesh0,
            const fvMesh& mesh1,
            const faceCoupleInfo& coupleInfo,
            const bool validBoundary = true,
            const bool fullyMapped = false
        );

        //- Map all volFields of Type
        template<class Type>
        static void MapVolFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd,
            const bool fullyMapped = false
        );

        //- Map all surfaceFields of Type
        template<class Type>
        static void MapSurfaceFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd,
            const bool fullyMapped = false
        );

        //- Map all surfaceFields of Type
        template<class Type>
        static void MapPointFields
        (
            const mapAddedPolyMesh&,
            const pointMesh& mesh,
            const labelListList& oldMeshPoints,
            const objectRegistry& meshToAdd
        );

        //- Map all DimensionedFields of Type
        template<class Type>
        static void MapDimFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshDistribute/meshAdderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

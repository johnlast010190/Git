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
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::parFvFieldReconstructor

Description
    Finite volume reconstructor for volume and surface fields.

    Runs in parallel. Reconstructs from procMesh to baseMesh. baseMesh
    is non-zero cells on processor0 only.

SourceFiles
    parFvFieldReconstructor.C
    parFvFieldReconstructorReconstructFields.C

\*---------------------------------------------------------------------------*/

#ifndef parFvFieldReconstructor_H
#define parFvFieldReconstructor_H

#include "containers/Lists/PtrList/PtrList.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyDistributionMap;
class distributionMapBase;
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class parFvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class parFvFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        fvMesh& baseMesh_;

        //- Processor mesh reference
        const fvMesh& procMesh_;

        //- Distribution map reference
        const polyDistributionMap& distMap_;

        //- Do I need to to write (usually master only)
        const bool isWriteProc_;

        //- Patch mappers
        PtrList<distributionMapBase> patchFaceMaps_;


    // Private Member Functions

        //- Construct per-patch addressing
        void createPatchFaceMaps();


public:

    // Constructors

        //- Construct from components
        parFvFieldReconstructor
        (
            fvMesh& baseMesh,
            const fvMesh& procMesh,
            const polyDistributionMap& distMap,
            const bool isWriteProc
        );

        //- Disallow default bitwise copy construction
        parFvFieldReconstructor(const parFvFieldReconstructor&) = delete;


    // Member Functions

        //- Reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        reconstructFvVolumeInternalField
        (
            const DimensionedField<Type, volMesh>&
        ) const;

        //- Read and reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        reconstructFvVolumeInternalField(const IOobject& fieldIoObject) const;



        //- Reconstruct volume field
        template<class Type>
        tmp<VolField<Type>>
        reconstructFvVolumeField
        (
            const VolField<Type>& fld
        ) const;

        //- Read and reconstruct volume field
        template<class Type>
        tmp<VolField<Type>>
        reconstructFvVolumeField(const IOobject& fieldIoObject) const;



        //- Reconstruct surface field
        template<class Type>
        tmp<SurfaceField<Type>>
        reconstructFvSurfaceField
        (
            const SurfaceField<Type>&
        ) const;

        //- Read and reconstruct surface field
        template<class Type>
        tmp<SurfaceField<Type>>
        reconstructFvSurfaceField(const IOobject& fieldIoObject) const;


        //- Read, reconstruct and write all/selected volume internal fields
        template<class Type>
        void reconstructFvVolumeInternalFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;

        //- Read, reconstruct and write all/selected volume fields
        template<class Type>
        void reconstructFvVolumeFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;

        //- Read, reconstruct and write all/selected surface fields
        template<class Type>
        void reconstructFvSurfaceFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;

        //- Helper: reconstruct and write mesh points
        //  (note: should be moved to something like processorMeshes class)
        void reconstructPoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const parFvFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "parFvFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

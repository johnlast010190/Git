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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::pointFieldReconstructor

Description
    Point field reconstructor.

SourceFiles
    pointFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldReconstructor_H
#define pointFieldReconstructor_H

#include "meshes/pointMesh/pointMesh.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"
#include "fields/Fields/fieldMappers/setSizeFieldMapper/setSizeFieldMapper.H"
#include "db/IOobjectList/IOobjectList.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class pointFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        const pointMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<pointMesh>& procMeshes_;

        //- List of processor point addressing lists
        const labelListList& pointProcAddressing_;

        //- Point patch addressing
        labelListListList patchPointAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


public:

        //- Mapper for sizing only - does not do any actual mapping.
        class pointPatchFieldReconstructor
        :
            public pointPatchFieldMapper,
            public setSizeFieldMapper
        {
        public:

            // Constructors

                //- Construct given size
                pointPatchFieldReconstructor(const label size)
                :
                    setSizeFieldMapper(size)
                {}
        };


    // Constructors

        //- Construct from components
        pointFieldReconstructor
        (
            const pointMesh& mesh,
            const PtrList<pointMesh>& procMeshes,
            const labelListList& pointProcAddressing
        );

        //- Disallow default bitwise copy construction
        pointFieldReconstructor(const pointFieldReconstructor&) = delete;


    // Member Functions

        //- Return number of fields reconstructed
        label nReconstructed() const
        {
            return nReconstructed_;
        }

        //- Reconstruct field
        template<class Type>
        tmp<PointField<Type>>
        reconstructField(const IOobject& fieldIoObject);

        //- Reconstruct and write all fields
        template<class Type>
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

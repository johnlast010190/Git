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
    (c) 2017 OpenCFD Ltd.

Class
    Foam::vtk::surfaceMeshWriter

Description
    Write faces (eg face-zones or face-sets) with fields.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkSurfaceMeshWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkSurfaceMeshWriter_H
#define foamVtkSurfaceMeshWriter_H

#include <fstream>
#include "fields/areaFields/areaFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "output/foamVtkIndPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::surfaceMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshWriter
:
    public vtk::indirectPatchWriter
{
    // Private Member Functions

        //- Get face field (internal face or boundary face)
        template<class Type>
        tmp<Field<Type>> getFaceField
        (
            const SurfaceField<Type>& sfld
        ) const;


        //- No copy construct
        surfaceMeshWriter(const surfaceMeshWriter&) = delete;

        //- No copy assignment
        void operator=(const surfaceMeshWriter&) = delete;


public:

    // Constructors

        //- Construct from patch (default format INLINE_BASE64)
        surfaceMeshWriter
        (
            const indirectPrimitivePatch& pp,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        )
        :
            vtk::indirectPatchWriter(pp, opts)
        {}

        //- Construct from patch (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        surfaceMeshWriter
        (
            const indirectPrimitivePatch& pp,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::indirectPatchWriter(pp, file, parallel)
        {}

        //- Construct from patch and open the file for writing.
        //  The file name is with/without an extension.
        surfaceMeshWriter
        (
            const indirectPrimitivePatch& pp,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::indirectPatchWriter(pp, opts, file, parallel)
        {}


    //- Destructor
    virtual ~surfaceMeshWriter() = default;


    // Member Functions

        //- Write surface field (CellData)
        template<class Type>
        void write
        (
            const SurfaceField<Type>& field
        );

        //- Write surface field (CellData)
        template<class Type>
        void write
        (
            const GeometricField<Type, faPatchField, areaMesh>& field
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtk/output/foamVtkSurfaceMeshWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

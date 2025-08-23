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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2017-2024 Engys Ltd.

Class
    Foam::sampledSurfaces

Description
    Set of surfaces to sample.

    The write() method is used to sample and write files.

    Example of function object specification:

    \verbatim
    surfaces
    {
        type    surfaces;
        libs    ("libsampling.so");

        // Write at same frequency as fields
        writeControl    outputTime;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Interpolation scheme to use (only used if interpolate=true for
        // the surfaces below)
        interpolationScheme cell;

        // Output surface format
        surfaceFormat   vtk;
        formatOptions   {  }

        surfaces
        (
            f0surf
            {
                type        sampledTriSurfaceMesh;
                surface     f0surf.obj;
                source      cells;

                // Optional: keep original regions
                keepIds     true;

                // Optional: generate values on points instead of faces
                interpolate true;
            }
        );
    }
    \endverbatim

See also
    Foam::surfMeshSamplers

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaces_H
#define sampledSurfaces_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "sampledSurface/writers/surfaceWriter.H"
#include "mergedSurf/mergedSurf.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/strings/lists/wordReList.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public functionObjects::regionFunctionObject,
    public PtrList<sampledSurface>
{
    // Static data members

        //- Output verbosity
        static bool verbose_;

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;


    // Private data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output path
        fileName outputPath_;

        // Read from dictonary

            //- Names of fields to sample
            wordReList fieldSelection_;

            //- Interpolation scheme to use
            word interpolationScheme_;

            //- whether to include time in output path
            bool writeTimeName_;

            //- Write field statistics
            bool writeStats_;

            //- Write sample surface fields
            bool writeFields_;

        // surfaces

            //- Merged meshed surfaces (parallel only)
            List<mergedSurf> mergedList_;


        // Calculated

            //- Surface formatter
            autoPtr<surfaceWriter> formatter_;


    // Private Member Functions


        //- Return number of fields
        label classifyFields();

        //- Write geometry only
        void writeGeometry() const;

        //- Write scalar field with original ids
        void writeOriginalIds();

         //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const VolField<Type>&
        );

        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const SurfaceField<Type>&
        );

        //- Sample and write all sampled fields
        template<class Type> void sampleAndWrite(const IOobjectList& objects);

         //- Sample and set a particular volume field
        template<class Type>
        void sampleAndSet
        (
            const VolField<Type>&,
            List<Field<Type>>& field
        );

        //- Sample and set a particular surface field
        template<class Type>
        void sampleAndSet
        (
            const SurfaceField<Type>&,
            List<Field<Type>>& field
        );

        //- Sample and set all sampled fields
        template<class GeoField, class Type>
        void sampleAndSet
        (
            const IOobjectList& objects,
            List<Field<Type>>& field
        );


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaces
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        sampledSurfaces(const sampledSurfaces&) = delete;


    //- Destructor
    virtual ~sampledSurfaces();


    // Member Functions

        //- Does any of the surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        //- Update the surfaces as required and merge surface points (parallel).
        //  Return false if no surfaces required an update.
        virtual bool update();


        //- Write sampled fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            const Field<Type>& values,
            const label surfI,
            const word& fieldName,
            const fileName& outputDir
        );

        //- Return whether to write fields
        virtual const bool& writeFields() const
        {
            return writeFields_;
        }

        //- Return whether to write field statistics
        virtual const bool& writeStats() const
        {
            return writeStats_;
        }
        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map - expires the surfaces
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update topology using the given map due to readUpdate
        //  - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set tolerance (and return old tolerance)
        static scalar mergeTol(const scalar);

        //- Sample and set field
        template<class GeoField, class Type>
        void returnField(List<Field<Type>>& field);

        //- Sample and set field
        template<class GeoField, class Type>
        Foam::List<Type> returnAverage(List<Field<Type>>& field);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sampledSurfaces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurface/sampledSurfaces/sampledSurfacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

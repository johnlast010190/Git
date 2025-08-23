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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::patchProbes

Description
    Set of locations to sample at patches

    Call write() to sample and write files.
    - find nearest location on nearest face
    - update *this with location (so header contains 'snapped' locations
    - use *this as the sampling location

    Example of function object specification:
    \verbatim
    patchProbes
    {
        type            patchProbes;
        libs            ( "libsampling.so" );

        // Name of the directory for probe data
        name            patchProbes;

        // Patches to sample (wildcards allowed)
        patches         (".*inl.*");

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be probed
        fields          (p U);

        // Locations to probe. These get snapped onto the nearest point
        // on the selected patches
        probeLocations
        (
            ( -100 0 0.01 )      // at inlet
        );
    }
    \endverbatim


SourceFiles
    patchProbes.C

\*---------------------------------------------------------------------------*/

#ifndef patchProbes_H
#define patchProbes_H

#include "probes/probes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class fvMesh;
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                          Class patchProbes Declaration
\*---------------------------------------------------------------------------*/

class patchProbes
:
    public probes
{
protected:

    // Protected data

        //- Patches to sample
        wordReList patchNames_;


    // Protected Member Functions

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

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

         //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const VolField<Type>&
        ) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const SurfaceField<Type>&
        ) const;

        //- Sample a single field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Find elements containing patchProbes
        virtual void findElements(const fvMesh&);


public:

    //- Runtime type information
    TypeName("patchProbes");


    // Constructors

        //- Construct from Time and dictionary
        patchProbes
        (
            const word& name,
            const Time& time,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

        //- Disallow default bitwise copy construction
        patchProbes(const patchProbes&) = delete;


    //- Destructor
    virtual ~patchProbes();


    //- Public members

        //- Execute, sample and write
        virtual bool execute();

        //- Write probe locations if coorFrame isIncrementalMotion
        virtual bool write();

        //- Read
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "probes/patchProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

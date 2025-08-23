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
    (c) 2011 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

Class
    Foam::multiFacePatchProbes

Description
    Use sphere to locate patch faces near probe point
    All faces contribute face area weighted properties to probe

    Call write() to sample and write files.

SourceFiles
    multiFacePatchProbes.C

\*---------------------------------------------------------------------------*/

#ifndef multiFacePatchProbes_H
#define multiFacePatchProbes_H

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
                          Class multiFacePatchProbes Declaration
\*---------------------------------------------------------------------------*/

class multiFacePatchProbes
:
    public probes
{
    // Private data

        //- Patch name
        word patchName_;

        //- probe size
        scalar radius2_;

        scalar probeScale_;

        //- face lists
        labelListList probeFaces_;

        //- probe surface sum for normalisation
        scalarField probeArea_;


    // Private Member Functions

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite(const VolField<Type>&);

        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite(const SurfaceField<Type>&);

        //- Sample and write all the fields of the given type
        template <class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

        //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const VolField<Type>&) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const SurfaceField<Type>&) const;

        //- Sample a single field on all sample locations
        template <class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;


public:

    //- Runtime type information
    TypeName("multiFacePatchProbes");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        multiFacePatchProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

        //- Disallow default bitwise copy construction
        multiFacePatchProbes(const multiFacePatchProbes&) = delete;


    //- Destructor
    virtual ~multiFacePatchProbes();


    // Public Member Functions

        //- Sample and write
        virtual bool write();

        //- Read
        virtual bool read(const dictionary&);

        //- Find elements containing multiFacePatchProbes
        virtual void findElements(const fvMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiFacePatchProbes&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "probes/multiFacePatchProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

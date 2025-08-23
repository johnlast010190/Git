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
    Foam::cloud

Description
    A cloud is a collection of lagrangian particles

SourceFiles
    cloud.C

\*---------------------------------------------------------------------------*/

#ifndef cloud_H
#define cloud_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/IOobjects/IOField/IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public objectRegistry
{
public:

    //- Runtime type information
    TypeName("cloud");

        //- The prefix to local: %lagrangian
        static const word prefix;

        //- The default cloud name: %defaultCloud
        static word defaultName;


    // Constructors

        //- Construct for the given objectRegistry and named cloud instance
        cloud(const objectRegistry&, const word& cloudName = defaultName);

        //- Disallow default bitwise copy construction
        cloud(const cloud&) = delete;


    //- Destructor
    virtual ~cloud();


    // Member Functions

        // Edit

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            virtual void autoMap(const polyTopoChangeMap&);


        // I-O

            //- Read particle fields from objects in the obr registry
            //virtual void readObjects(objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            virtual void writeObjects(objectRegistry& obr) const;

            //- Helper to construct IOField on a supplied object registry
            template<class Type>
            static IOField<Type>& createIOField
            (
                const word& fieldName,
                const label nParticle,
                objectRegistry& obr
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/cloud/cloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

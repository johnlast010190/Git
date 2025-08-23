/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::IOOutputFilter

Description
    IOdictionary wrapper around OutputFilter to allow them to read from
    their associated dictionaries.

Note
    The IOobject or the objectRegistry will normally have to be
    derived from a fvMesh for a subsequent cast (within OutputFilter)
    to work correctly.

SourceFiles
    IOOutputFilter.C

\*---------------------------------------------------------------------------*/

#ifndef IOOutputFilter_H
#define IOOutputFilter_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "meshes/primitiveShapes/point/pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;
class polyMeshMap;

/*---------------------------------------------------------------------------*\
                       Class IOOutputFilter Declaration
\*---------------------------------------------------------------------------*/

template<class OutputFilter>
class IOOutputFilter
:
    public IOdictionary,
    public OutputFilter
{
public:

    // Constructors

        //- Construct from an IOobject for IOdictionary
        //  Allow the possibility to load fields from files
        IOOutputFilter
        (
            const word& outputFilterName,
            const IOobject& ioDict,
            const bool loadFromFile = false
        );

        //- Construct for given objectRegistry and dictionary
        //  Allow dictionary to be optional
        //  Allow the possibility to load fields from files
        IOOutputFilter
        (
            const word& outputFilterName,
            const objectRegistry&,
            const word& dictName = OutputFilter::typeName() + "Dict",
            const IOobject::readOption rOpt = IOobject::MUST_READ_IF_MODIFIED,
            const bool loadFromFile = false
        );

        //- Construct for given objectRegistry and dictionary
        //  Dictionary read from full path.
        //  Allow the possibility to load fields from files
        IOOutputFilter
        (
            const word& outputFilterName,
            const objectRegistry&,
            const fileName& dictName,
            const IOobject::readOption rOpt = IOobject::MUST_READ_IF_MODIFIED,
            const bool loadFromFile = false
        );

        //- Disallow default bitwise copy construction
        IOOutputFilter(const IOOutputFilter&) = delete;


    //- Destructor
    virtual ~IOOutputFilter();


    // Member Functions

        //- Return name
        virtual const word& name() const
        {
            return IOdictionary::name();
        }

        //- Inherit read from OutputFilter
        using OutputFilter::read;

        //- Read output filter properties
        virtual bool read();

        //- Inherit write from regIOobject
        using regIOobject::write;

        //- Sample and write
        virtual bool write();

        //- Update topology using the given map
        virtual void movePoints(const polyMesh& mesh);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IOOutputFilter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOOutputFilter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::IOMap

Description
    A Map of objects of type \<T\> with automated input and output.
    Is a global object; i.e. can be read from undecomposed case.

SourceFiles
    IOMap.C

\*---------------------------------------------------------------------------*/

#ifndef IOMap_H
#define IOMap_H

#include "containers/HashTables/Map/Map.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOMap
:
    public regIOobject,
    public Map<T>
{

public:

    //- Runtime type information
    TypeName("Map");


    // Constructors

        //- Construct from IOobject
        IOMap(const IOobject&);

        //- Construct from IOobject and size of Map
        IOMap(const IOobject&, const label);

        //- Construct from IOobject and a Map
        IOMap(const IOobject&, const Map<T>&);

        //- Move construct by transferring the Map contents
        IOMap(const IOobject&, Map<T>&&);

        //- Move constructor
        IOMap(IOMap<T>&&);


    //- Destructor
    virtual ~IOMap();


    // Member functions

        bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }


    // Member operators

        //- Assignment of other IOMap's entries to this IOMap
        void operator=(const IOMap<T>&);

        //- Move assignment of other IOMap's entries to this IOMap
        void operator=(IOMap<T>&&);

        //- Assignment of other Map's entries to this IOMap
        void operator=(const Map<T>&);

        //- Move assignment of other Map's entries to this IOMap
        void operator=(Map<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/IOobjects/IOMap/IOMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

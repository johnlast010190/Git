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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::IOdictionary

Description
    IOdictionary is derived from dictionary and IOobject to give the dictionary
    automatic IO functionality via the objectRegistry.

    To facilitate IO, IOdictionary is provided with a constructor from IOobject
    and with readData/writeData functions.

SourceFiles
    IOdictionary.C
    IOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef IOdictionary_H
#define IOdictionary_H

#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IOdictionary Declaration
\*---------------------------------------------------------------------------*/

class IOdictionary
:
    public regIOobject,
    public dictionary
{
    // Private Data

        static bool writeDictionaries;


protected:

    // Protected constructors

        //- Construct given an IOobject, supply wanted typeName
        IOdictionary(const IOobject& io, const word& wantedType);


public:

    TypeName("dictionary");


    // Constructors

        //- Construct given an IOobject
        IOdictionary(const IOobject&);

        //- Construct given an IOobject and dictionary
        IOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and Istream
        IOdictionary(const IOobject&, Istream&);

        //- Copy constructor
        IOdictionary(const IOdictionary&);

        //- Move constructor
        IOdictionary(IOdictionary&&);


    //- Destructor
    virtual ~IOdictionary();


    // Member functions

        using regIOobject::name;

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

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


// Member operators

        //- Assignment of dictionary entries only
        void operator=(const IOdictionary&);

        //- Move assignment
        void operator=(IOdictionary&&);
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<IOdictionary>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

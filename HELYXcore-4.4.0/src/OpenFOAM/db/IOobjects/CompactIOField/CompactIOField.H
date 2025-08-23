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
    Foam::CompactIOField

Description
    A Field of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOField except when binary output in
    case it writes a CompactListList.

    Useful for fields of small subfields e.g. in lagrangian

SourceFiles
    CompactIOField.C

\*---------------------------------------------------------------------------*/

#ifndef CompactIOField_H
#define CompactIOField_H

#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T, class BaseType> class CompactIOField;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOField<T, BaseType>&
);

template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOField<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                        Class CompactIOField Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOField
:
    public regIOobject,
    public Field<T>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream(const bool read = true);

public:

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Construct from IOobject
        CompactIOField(const IOobject&);

        //- Construct from IOobject; does local processor require reading?
        CompactIOField(const IOobject&, const bool read);

        //- Construct from IOobject and size
        CompactIOField(const IOobject&, const label);

        //- Construct from IOobject and a Field
        CompactIOField(const IOobject&, const Field<T>&);

        //- Move construct by transferring the Field contents
        CompactIOField(const IOobject&, Field<T>&&);

        //- Move constructor
        CompactIOField(CompactIOField<T, BaseType>&&);


    // Destructor

        virtual ~CompactIOField();


    // Member functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member operators

        void operator=(const CompactIOField<T, BaseType>&);
        void operator=(CompactIOField<T, BaseType>&&);

        void operator=(const Field<T>&);
        void operator=(Field<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/IOobjects/CompactIOField/CompactIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

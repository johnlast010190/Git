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
    (c) 2017 OpenCFD Ltd.
    (c) 2014-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::nullObject

Description
    Singleton null-object class and instance.
    It occupies enough space to reinterpret its content as a class with
    a null pointer for its content.

SourceFiles
    nullObjectI.H
    nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef nullObject_H
#define nullObject_H

#include <utility>  // for std::move

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class nullObject Declaration
\*---------------------------------------------------------------------------*/

class NullObject
{
    //- Ensure it occupies enough space to reinterpret_cast to a class
    //  having some member data
    const union
    {
        void* ptr;
        unsigned long val;
    } null;

    //- Private constructor
    NullObject()
    :
        null{nullptr}
    {}

    //- Disallow default bitwise copy construction
    NullObject(const NullObject&) = delete;


public:

    //- The unique null object
    static const NullObject nullObject;

    //- A nullptr pointer content
    inline const void* pointer() const
    {
        return null.ptr;
    }

    //- A zero value content
    inline unsigned long value() const
    {
        return null.val;
    }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const NullObject&) = delete;
};


//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


//- Return const reference to the nullObject of type T
template<class T>
inline const T& NullObjectRef();

//- Return const pointer to the nullObject of type T
template<class T>
inline const T* NullObjectConstPtr();

//- Return pointer to the nullObject of type T
template<class T>
inline T* NullObjectPtr();

//- Return reference to the nullObject of type T
template<class T>
inline T&& NullObjectMove();

//- Return true if t is a reference to the nullObject of type T
template<class T>
inline bool isNull(const T& t);

//- Return true if t is not a reference to the nullObject of type T
template<class T>
inline bool notNull(const T& t);


//- Return true if t is a pointer to the nullObject of type T
template<class T>
inline bool isNull(const T* t);

//- Return true if t is not a pointer to the nullObject of type T
template<class T>
inline bool notNull(const T* t);


//- Return true if t is a pointer to the nullObject of type T
template<class T>
inline bool isNull(T* t);

//- Return true if t is not a pointer to the nullObject of type T
template<class T>
inline bool notNull(T* t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/nullObject/nullObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

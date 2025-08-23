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
    Foam::PtrDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on a memory managed form of
    intrusive doubly-linked list of \<T\>.

SourceFiles
    PtrDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef PtrDictionary_H
#define PtrDictionary_H

#include "containers/Dictionaries/DictionaryBase/DictionaryBase.H"
#include "containers/LinkedLists/user/DLPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrDictionary
:
    public DictionaryBase<DLPtrList<T>, T>
{

public:

    // Constructors

        //- Construct given initial table size
        PtrDictionary(const label size = 128);

        //- Copy construct
        PtrDictionary(const PtrDictionary&);

        //- Move constructor
        PtrDictionary(PtrDictionary&&);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrDictionary(Istream&, const INew&);

        //- Construct from Istream
        PtrDictionary(Istream&);


    // Member operators

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *DictionaryBase<DLPtrList<T>, T>::operator[](key);
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *DictionaryBase<DLPtrList<T>, T>::operator[](key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "containers/Dictionaries/PtrDictionary/PtrDictionary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

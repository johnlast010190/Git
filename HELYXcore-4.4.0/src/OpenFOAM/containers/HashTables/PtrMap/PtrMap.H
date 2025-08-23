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
    (c) 2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::PtrMap

Description
    A HashTable of pointers to objects of type \<T\> with a label key.

See also
    Map, HashPtrTable

\*---------------------------------------------------------------------------*/

#ifndef PtrMap_H
#define PtrMap_H

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PtrMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrMap
:
    public HashPtrTable<T, label, Hash<label>>
{
public:

    //- The template instance used for this PtrMap
    typedef PtrMap<T> this_type;

    //- The template instance used for the parent HashTable
    typedef HashPtrTable<T, label, Hash<label>> parent_type;


    // Constructors

        //- Construct given initial map size
        PtrMap(const label size = 128)
        :
            parent_type(size)
        {}

        //- Construct from Istream
        PtrMap(Istream& is)
        :
            parent_type(is)
        {}

        //- Copy constructor
        PtrMap(const this_type& map)
        :
            parent_type(map)
        {}

        //- Move constructor
        PtrMap(PtrMap<T>&& map)
        :
            parent_type(std::move(map))
        {}


    // Member operators

        //- Copy assignment operator
        void operator=(const PtrMap<T>& map)
        {
            parent_type::operator=(map);
        }

        //- Move assignment operator
        void operator=(PtrMap<T>&& map)
        {
            parent_type::operator=(std::move(map));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

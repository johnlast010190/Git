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
    Foam::Map

Description
    A HashTable to objects of type \<T\> with a label key.

See also
    PtrMap

\*---------------------------------------------------------------------------*/

#ifndef Map_H
#define Map_H

#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Map Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Map
:
    public HashTable<T, label, Hash<label>>
{
public:

    //- The template instance used for this Map
    typedef Map<T> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<T, label, Hash<label>> parent_type;

    using iterator = typename parent_type::iterator;
    using const_iterator = typename parent_type::const_iterator;


    // Constructors

        //- Construct given initial size
        Map(const label size = 128)
        :
            parent_type(size)
        {}

        //- Construct from Istream
        Map(Istream& is)
        :
            parent_type(is)
        {}

        //- Copy constructor
        Map(const Map<T>& map)
        :
            parent_type(map)
        {}

        //- Move constructor
        Map(Map<T>&& map)
        :
            parent_type(std::move(map))
        {}

        //- Move constructor
        Map(parent_type&& map)
        :
            parent_type(std::move(map))
        {}

        //- Construct from an initializer list
        Map(std::initializer_list<std::pair<label, T>> map)
        :
            parent_type(map)
        {}

    // Member Operators

        void operator=(const Map<T>& map)
        {
            parent_type::operator=(map);
        }

        void operator=(Map<T>&& map)
        {
            parent_type::operator=(std::move(map));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

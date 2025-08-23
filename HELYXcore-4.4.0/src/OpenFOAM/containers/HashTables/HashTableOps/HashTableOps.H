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
    (c) 2018 OpenCFD Ltd.

InNamspace
    Foam

Description
    Various functions to operate on HashTables.

SourceFiles
    HashTableOps.H

\*---------------------------------------------------------------------------*/

#ifndef HashTableOps_H
#define HashTableOps_H

#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Combine HashSet operation. Equivalent to 'a += b'
template<class Key=word, class Hash=string::hash>
struct HashSetPlusEqOp
{
    typedef HashSet<Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};


//- Combine HashTable operation. Equivalent to 'a += b'
template<class T, class Key=word, class Hash=string::hash>
struct HashTablePlusEqOp
{
    typedef HashTable<T, Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        if (b.size())
        {
            if (a.size())
            {
                forAllConstIters(b, citer)
                {
                    a.insert(citer.key(), citer.object());
                }
            }
            else
            {
                a = b;
            }
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

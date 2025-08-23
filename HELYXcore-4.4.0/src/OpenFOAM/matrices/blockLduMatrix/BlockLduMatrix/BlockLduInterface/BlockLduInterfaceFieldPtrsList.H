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
    (c) 2004-2006 Ivor Clifford

Class
    BlockLduInterfaceFieldPtrsList

Description
    List of coupled interface fields to be used in coupling.
    Templated typedef workaround.

Author
    Ivor Clifford, ivor.clifford@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterfaceFieldPtrsList_H
#define BlockLduInterfaceFieldPtrsList_H

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduInterface/BlockLduInterfaceField.H"
#include "containers/Lists/UPtrList/UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
struct BlockLduInterfaceFieldPtrsList
{
    typedef UPtrList<const BlockLduInterfaceField<T>> Type;
};

template<class T>
struct BlockLduInterfaceFieldPtrsListList
{
    typedef List<typename BlockLduInterfaceFieldPtrsList<T>::Type> Type;
};

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "containers/LinkedLists/accessTypes/LPtrList/LPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(const LPtrList<LListBase, T>& lst)
{
    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(iter().clone().ptr());
    }
}


template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(LPtrList<LListBase, T>&& lst)
{
    transfer(lst);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::~LPtrList()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LListBase, class T>
bool Foam::LPtrList<LListBase, T>::eraseHead()
{
    T* tPtr;
    if ((tPtr = this->removeHead()))
    {
        delete tPtr;
        return true;
    }
    else
    {
        return false;
    }
}


template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::clear()
{
    const label oldSize = this->size();
    for (label i=0; i<oldSize; ++i)
    {
        eraseHead();
    }

    LList<LListBase, T*>::clear();
}


template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::transfer(LPtrList<LListBase, T>& lst)
{
    clear();
    LList<LListBase, T*>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::operator=(const LPtrList<LListBase, T>& lst)
{
    clear();

    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(iter().clone().ptr());
    }
}


template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::operator=(LPtrList<LListBase, T>&& lst)
{
    transfer(lst);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "containers/LinkedLists/accessTypes/LPtrList/LPtrListIO.C"


// ************************************************************************* //

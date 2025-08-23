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

    (c) 2018-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::IOobject* Foam::IOobjectList::cfindObject
(
    const word& objName
) const
{
    const_iterator iter = cfind(objName);

    if (iter.found())
    {
        const IOobject* io = iter();

        if (io->isHeaderClassName<Type>())
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << objName << endl;
            }

            return io;
        }
        else if (IOobject::debug)
        {
            InfoInFunction
                << "Found " << objName << " of different type" << endl;
        }
    }
    else if (IOobject::debug)
    {
        InfoInFunction << "Could not find " << objName << endl;
    }

    return nullptr;
}


// Templated implementation for names(), sortedNames()
template<class Type>
Foam::wordList Foam::IOobjectList::namesTypeImpl
(
    const IOobjectList& list,
    const bool doSort
)
{
    wordList objNames(list.size());

    label count = 0;
    forAllConstIters(list, iter)
    {
        const word& key = iter.key();
        const IOobject* io = iter();

        if (io->isHeaderClassName<Type>())
        {
            objNames[count] = key;
            ++count;
        }
    }

    objNames.resize(count);

    if (doSort)
    {
        Foam::sort(objNames);
    }

    return objNames;
}



template<class Type>
Foam::wordList Foam::IOobjectList::sortedNames() const
{
    // return namesTypeImpl<Type>(*this, predicates::always(), true);
    return namesTypeImpl<Type>(*this, true);
}

template<class Type>
Foam::wordList Foam::IOobjectList::allNames() const
{
    // wordList objNames(namesTypeImpl<Type>(*this, predicates::always(), false));
    wordList objNames(namesTypeImpl<Type>(*this, false));

    syncNames(objNames);
    return objNames;
}


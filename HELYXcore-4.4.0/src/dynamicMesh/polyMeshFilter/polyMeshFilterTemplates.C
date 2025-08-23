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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "polyMeshFilter/polyMeshFilter.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class SetType>
void Foam::polyMeshFilter::updateSets(const polyTopoChangeMap& map)
{
    HashTable<const SetType*> sets =
        map.mesh().objectRegistry::lookupClass<const SetType>();

    forAllIter(typename HashTable<const SetType*>, sets, iter)
    {
        SetType& set = const_cast<SetType&>(*iter());
        set.topoChange(map);
        set.sync(map.mesh());
    }

    IOobjectList Objects
    (
        map.mesh().time(),
        map.mesh().facesInstance(),
        "polyMesh/sets"
    );

    IOobjectList fileSets(Objects.lookupClass(SetType::typeName));

    forAllConstIter(IOobjectList, fileSets, iter)
    {
        if (!sets.found(iter.key()))
        {
            // Not in memory. Load it.
            SetType set(*iter());
            set.topoChange(map);

            set.write();
        }
    }
}


template<class SetType>
void Foam::polyMeshFilter::copySets
(
    const polyMesh& oldMesh,
    const polyMesh& newMesh
)
{
    HashTable<const SetType*> sets =
        oldMesh.objectRegistry::lookupClass<const SetType>();

    forAllConstIter(typename HashTable<const SetType*>, sets, iter)
    {
        const SetType& set = *iter();

        if (newMesh.objectRegistry::foundObject<SetType>(set.name()))
        {
            const SetType& origSet =
                newMesh.objectRegistry::lookupObject<SetType>(set.name());

            const_cast<SetType&>(origSet) = set;
            const_cast<SetType&>(origSet).sync(newMesh);
        }
        else
        {
            SetType* newSet
            (
                new SetType(newMesh, set.name(), set, set.writeOpt())
            );

            newSet->store();
            newSet->sync(newMesh);
        }
    }
}


// ************************************************************************* //

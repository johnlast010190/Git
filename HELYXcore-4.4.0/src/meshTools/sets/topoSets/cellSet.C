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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyDistributionMap/polyDistributionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSet, 0);

addToRunTimeSelectionTable(topoSet, cellSet, word);
addToRunTimeSelectionTable(topoSet, cellSet, size);
addToRunTimeSelectionTable(topoSet, cellSet, set);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::cellSet::cellSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    // Make sure set within valid range
    check(mesh.nCells());
}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const UList<label>& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


// Database constructors (for when no mesh available)
Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, r, w),
        typeName
    )
{}


Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, IOobject::NO_READ, w),
        size
    )
{}


Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, IOobject::NO_READ, w),
        set
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSet::~cellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellSet::topoChange(const polyTopoChangeMap& map)
{
    updateLabels(map.reverseCellMap());
}


void Foam::cellSet::distribute(const polyDistributionMap& map)
{
    boolList inSet(map.nOldCells(),false);
    forAllConstIter(cellSet, *this, iter)
    {
        inSet[iter.key()] = true;
    }
    map.distributeCellData(inSet);

    // Count
    label n = 0;
    forAll(inSet, celli)
    {
        if (inSet[celli])
        {
            n++;
        }
    }

    clear();
    resize(n);
    forAll(inSet, celli)
    {
        if (inSet[celli])
        {
            insert(celli);
        }
    }
}


void Foam::cellSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// ************************************************************************* //

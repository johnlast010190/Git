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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFvMesh, 0);
    defineRunTimeSelectionTable(dynamicFvMesh, IOobject);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::dynamicFvMesh::dynamicMeshDictIOobject(const IOobject& io)
{
    // defaultRegion (region0) gets loaded from constant, other ones get loaded
    // from constant/<regionname>. Normally we'd use polyMesh::dbDir() but we
    // haven't got a polyMesh yet...
    return IOobject
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::dynamicFvMesh(const IOobject& io)
:
    fvMesh(io, false, fvMesh::stitchType::none),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{
    DeprecationWarningInFunction
    (
        "dynamicFvMesh",
        "dynamic mesh library",
        40300,
        "Please replace it with the new multi-motion dynamic library, "
        "by setting the 'meshChangers' keyword to 'true' in 'controlDict' "
        "and using mesh changers instead."
    );
}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    ),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{}


bool Foam::dynamicFvMesh::update(const word& cfg)
{
    WarningInFunction
        << "Function not implemented for the " << this->typeName
        << " object type." << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::~dynamicFvMesh()
{}


// ************************************************************************* //

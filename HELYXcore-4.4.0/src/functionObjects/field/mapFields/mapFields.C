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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2018 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "mapFields/mapFields.H"
#include "meshToMesh/meshToMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mapFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mapFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mapFields::createInterpolation
(
    const dictionary& dict
)
{
    const fvMesh& meshTarget = mesh_;
    const word mapRegionName(dict.lookup("mapRegion"));

    Info<< name() << ':' << nl
        << "    Reading mesh " << mapRegionName << endl;

    mapRegionPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                mapRegionName,
                meshTarget.time().constant(),
                meshTarget.time()
            )
        )
    );
    const fvMesh& mapRegion = mapRegionPtr_();
    const word mapMethodName(dict.lookup("mapMethod"));
    if (!meshToMesh::interpolationMethodNames_.hasEnum(mapMethodName))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Lookup corresponding AMi method
    word patchMapMethodName =
        AMIInterpolation::interpolationMethodToWord
        (
            meshToMesh::interpolationMethodAMI(mapMethod)
        );

    // Optionally override
    if (dict.readIfPresent("patchMapMethod", patchMapMethodName))
    {
        Info<< "    Patch mapping method: " << patchMapMethodName << endl;
    }

    bool consistent = dict.lookup<bool>("consistent");

    Info<< "    Creating mesh to mesh interpolation" << endl;

    if (consistent)
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName
            )
        );
    }
    else
    {
        HashTable<word> patchMap(dict.lookup("patchMap"));
        wordList cuttingPatches(dict.lookup("cuttingPatches"));

        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName,
                patchMap,
                cuttingPatches
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::mapFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mapRegionPtr_(),
    interpPtr_(),
    fieldNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::~mapFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mapFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    fieldNames_ = dict.lookup<wordReList>("fields");
    createInterpolation(dict);
    return true;
}


bool Foam::functionObjects::mapFields::execute()
{
    return true;
}


bool Foam::functionObjects::mapFields::write()
{
    Log << type() << " " << name() << " write:" << nl;

    bool ok = false;

    ok = writeFieldType<scalar>() || ok;
    ok = writeFieldType<vector>() || ok;
    ok = writeFieldType<sphericalTensor>() || ok;
    ok = writeFieldType<symmTensor>() || ok;
    ok = writeFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }

    return true;
}


// ************************************************************************* //

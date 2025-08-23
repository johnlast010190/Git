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
    (c) 2022-2023 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshToMesh/fvMeshTopoChangersMeshToMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "fields/volFields/volFields.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "meshToMesh/meshToMeshAdjustTimeStep/meshToMeshAdjustTimeStepFunctionObject.H"
#include "meshToMesh/meshToMesh.H"
#include "meshToMesh/calcMethod/cellVolumeWeight/cellVolumeWeightMethod.H"
#include "cfdTools/general/surfaceToVolVelocity/surfaceToVolVelocity.H"
#include "meshToMesh/MeshToMeshMapGeometricFields.H"
#include "meshes/polyMesh/polyMeshMap/polyMeshMap.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(meshToMesh, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, meshToMesh, fvMesh);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::meshToMesh::interpolateUfs()
{
    // Interpolate U to Uf
    UPtrList<surfaceVectorField> Ufs(mesh().curFields<surfaceVectorField>());

    forAll(Ufs, i)
    {
        surfaceVectorField& Uf = Ufs[i];

        const volVectorField& U = surfaceToVolVelocity(Uf);

        if (!isNull(U))
        {
            Uf.reset(fvc::interpolate(U));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::meshToMesh::meshToMesh
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    times_(dict.lookup("times")),
    timeDelta_(dict.lookup<scalar>("timeDelta")),
    timeIndex_(-1)
{
    forAll(times_, i)
    {
        timeIndices_.insert(label(times_[i]/timeDelta_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::meshToMesh::update()
{
    if (timeIndex_ == -1)
    {
        const_cast<Time&>(mesh().time()).functionObjects().append
        (
            new functionObjects::meshToMeshAdjustTimeStepFunctionObject
            (
                "meshToMeshAdjustTimeStep",
                mesh().time(),
                dict_
            )
        );
    }

    bool hasChanged = false;

    // Only refine on the first call in a time-step
    if (timeIndex_ != mesh().time().timeIndex())
    {
        timeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        return hasChanged;
    }

    const scalar meshTime =
        mesh().time().timeToUserTime(mesh().time().value())
      - mesh().time().timeToUserTime(mesh().time().deltaTValue());

    if (timeIndices_.found((meshTime + timeDelta_/2)/timeDelta_))
    {
        const word otherMeshDir =
            "meshToMesh_" + mesh().time().timeName(meshTime);

        Info<< "Mapping to mesh " << otherMeshDir << endl;

        hasChanged = true;

        fvMesh otherMesh
        (
            IOobject
            (
                otherMeshDir,
                mesh().time().constant(),
                mesh().time(),
                IOobject::MUST_READ
            ),
            false,
            fvMesh::stitchType::none
        );

        mesh().swap(otherMesh);

        Foam::meshToMesh mapper
        (
            otherMesh,
            mesh(),
            Foam::meshToMesh::imCellVolumeWeight
        );

        // Ensure the deltaCoeffs are available for constraint patch evaluation
        mesh().deltaCoeffs();

        // Map all the volFields in the objectRegistry
        #define mapVolFieldType(Type, nullArg)                                 \
            MeshToMeshMapVolFields<Type>(mesh(), mapper);
        FOR_ALL_FIELD_TYPES(mapVolFieldType);
/*
        // Set all the surfaceFields in the objectRegistry to NaN
        #define NaNSurfaceFieldType(Type, nullArg)                             \
            NaNGeometricFields                                                 \
            <Type, fvsPatchField, surfaceMesh, fvPatchFieldMapper>             \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNSurfaceFieldType);
*/
        // Set all the pointFields in the objectRegistry to NaN
        #define NaNPointFieldType(Type, nullArg)                               \
            NaNGeometricFields                                                 \
            <Type, pointPatchField, pointMesh, pointPatchFieldMapper>          \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNPointFieldType);

        // Interpolate U's to Uf's
        interpolateUfs();

        polyMeshMap map(mesh());
        mesh().mapMesh(map);
    }

    return hasChanged;
}


void Foam::fvMeshTopoChangers::meshToMesh::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fvMeshTopoChangers::meshToMesh::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::meshToMesh::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //

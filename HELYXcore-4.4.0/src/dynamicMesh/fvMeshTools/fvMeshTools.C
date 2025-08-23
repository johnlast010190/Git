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
    (c) 2015 OpenCFD Ltd.
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshTools/fvMeshTools.H"
#include "meshes/polyMesh/polyPatches/constraint/processorCyclic/processorCyclicPolyPatch.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMeshEntries.H"
#include "fvMeshTools/fvMeshTools.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "include/OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Adds patch if not yet there. Returns patchID.
Foam::label Foam::fvMeshTools::addPatch
(
    fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const bool validBoundary
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchi = polyPatches.findPatchID(patch.name());
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    // Append at end unless there are processor patches
    label insertPatchi = polyPatches.size();

    if (!isA<processorPolyPatch>(patch))
    {
        forAll(polyPatches, patchi)
        {
            const polyPatch& pp = polyPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                insertPatchi = patchi;
                break;
            }
        }
    }

    mesh.addPatch
    (
        insertPatchi,
        patch,
        patchFieldDict,
        defaultPatchFieldType,
        validBoundary
    );

    return insertPatchi;
}


Foam::label Foam::fvMeshTools::addPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchDict,
    const bool validBoundary
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    // create patch from dict
    autoPtr<polyPatch> ppPtr
    (
        polyPatch::New
        (
            patchName,
            patchDict,
            polyPatches.size(),
            polyPatches
        )
    );

    // call addPatch
    return addPatch
    (
        mesh,
        ppPtr(),
        dictionary(), // maybe allow optional patchFieldDict
        fvPatchField<scalar>::calculatedType(),
        validBoundary
    );
}


void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    setPatchFields<volScalarField>(mesh, patchi, patchFieldDict);
    setPatchFields<volVectorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volSphericalTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volSymmTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volTensorField>(mesh, patchi, patchFieldDict);

    setPatchFields<surfaceScalarField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceVectorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceSphericalTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceSymmTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceTensorField>(mesh, patchi, patchFieldDict);

    setPatchFields<volVector1Field>(mesh, patchi, patchFieldDict);
    setPatchFields<volVector4Field>(mesh, patchi, patchFieldDict);
    setPatchFields<volTensor4Field>(mesh, patchi, patchFieldDict);

    if (mesh.foundObject<pointMesh>(pointMesh::typeName))
    {
        pointMesh& pm = const_cast<pointMesh&>(pointMesh::New(mesh));
        setPatchFields<pointScalarField>(pm, patchi, patchFieldDict);
        setPatchFields<pointVectorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointSphericalTensorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointSymmTensorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointTensorField>(pm, patchi, patchFieldDict);
    }
}


void Foam::fvMeshTools::zeroPatchFields(fvMesh& mesh, const label patchi)
{
    setPatchFields<volScalarField>(mesh, patchi, Zero);
    setPatchFields<volVectorField>(mesh, patchi, Zero);
    setPatchFields<volSphericalTensorField>(mesh, patchi, Zero);
    setPatchFields<volSymmTensorField>(mesh, patchi, Zero);
    setPatchFields<volTensorField>(mesh, patchi, Zero);

    setPatchFields<volVector1Field>(mesh, patchi, Zero);
    setPatchFields<volVector4Field>(mesh, patchi, Zero);
    setPatchFields<volTensor4Field>(mesh, patchi, Zero);

    setPatchFields<surfaceScalarField>(mesh, patchi, Zero);
    setPatchFields<surfaceVectorField>(mesh, patchi, Zero);
    setPatchFields<surfaceSphericalTensorField>(mesh, patchi, Zero);
    setPatchFields<surfaceSymmTensorField>(mesh,patchi,Zero);
    setPatchFields<surfaceTensorField>(mesh, patchi, Zero);

    if (mesh.foundObject<pointMesh>(pointMesh::typeName))
    {
        pointMesh& pm = const_cast<pointMesh&>(pointMesh::New(mesh));
        setPatchFields<pointScalarField>(pm, patchi, Zero);
        setPatchFields<pointVectorField>(pm, patchi, Zero);
        setPatchFields<pointSphericalTensorField>(pm, patchi, Zero);
        setPatchFields<pointSymmTensorField>(pm, patchi, Zero);
        setPatchFields<pointTensorField>(pm, patchi, Zero);
    }
}


void Foam::fvMeshTools::reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches,
    const bool validBoundary
)
{
    labelList newToOld(nNewPatches, -1);
    forAll(oldToNew, i)
    {
        label newi = oldToNew[i];

        if (newi >= 0 && newi < nNewPatches)
        {
            newToOld[newi] = i;
        }
    }

    mesh.reorderPatches(newToOld, validBoundary);
}


Foam::labelList Foam::fvMeshTools::removeEmptyPatches
(
    fvMesh& mesh,
    const bool validBoundary
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList newToOld(pbm.size());
    labelList oldToNew(pbm.size(), -1);
    label newI = 0;


    // Assumes all non-coupled boundaries are on all processors!
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (!isA<processorPolyPatch>(pp))
        {
            label nFaces = pp.size();
            if (validBoundary)
            {
                reduce(nFaces, sumOp<label>());
            }

            if (nFaces > 0)
            {
                newToOld[newI] = patchI;
                oldToNew[patchI] = newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            newToOld[newI] = patchI;
            oldToNew[patchI] = newI++;
        }
    }

    newToOld.setSize(newI);

    // Move all deleteable patches to the end
    forAll(oldToNew, patchI)
    {
        if (oldToNew[patchI] == -1)
        {
            oldToNew[patchI] = newI++;
        }
    }

    reorderPatches(mesh, oldToNew, newToOld.size(), validBoundary);

    return newToOld;
}


Foam::autoPtr<Foam::fvMesh> Foam::fvMeshTools::newMesh
(
    const IOobject& io,
    const bool masterOnlyReading
)
{
    // Region name
    // ~~~~~~~~~~~

    fileName meshSubDir;

    if (io.name() == polyMesh::defaultRegion)
    {
        meshSubDir = polyMesh::meshSubDir;
    }
    else
    {
        meshSubDir = io.name()/polyMesh::meshSubDir;
    }


    fileName facesInstance;

    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (Pstream::master())
    {
        facesInstance = io.time().findInstance
        (
            meshSubDir,
            "faces",
            IOobject::MUST_READ
        );

        patchEntries = polyBoundaryMeshEntries
        (
            IOobject
            (
                "boundary",
                facesInstance,
                meshSubDir,
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Send patches
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(Pstream::commsTypes::scheduled, slave);
            toSlave << patchEntries;
        }
    }
    else
    {
        // Receive patches
        IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
        fromMaster >> patchEntries;
    }


    Pstream::scatter(facesInstance);

    // Dummy meshes
    // ~~~~~~~~~~~~

    // Check who has a mesh
    const fileName meshDir = io.time().path()/facesInstance/meshSubDir;
    bool haveMesh = isDir(meshDir);
    if (masterOnlyReading && !Pstream::master())
    {
        haveMesh = false;
    }


    // Set up to read-if-present. Note: does not search for mesh so set
    // instance explicitly
    IOobject meshIO(io);
    meshIO.instance() = facesInstance;
    if (masterOnlyReading && !Pstream::master())
    {
        meshIO.readOpt() = IOobject::NO_READ;
    }
    else
    {
        meshIO.readOpt() = IOobject::READ_IF_PRESENT;
    }


    // Read mesh
    // ~~~~~~~~~
    // Now all processors read a mesh and use supplied points,faces etc
    // if there is none.
    // Note: fvSolution, fvSchemes are also using the supplied Ioobject so
    //       on slave will be NO_READ, on master READ_IF_PRESENT. This will
    //       conflict with e.g. timeStampMaster reading so switch off.

    const regIOobject::fileCheckTypes oldCheckType =
        regIOobject::fileModificationChecking;
    regIOobject::fileModificationChecking = regIOobject::timeStamp;

    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            meshIO,
            pointField(),
            faceList(),
            labelList(),
            labelList()
        )
    );
    fvMesh& mesh = meshPtr();

    regIOobject::fileModificationChecking = oldCheckType;



    // Add patches
    // ~~~~~~~~~~~


    bool havePatches = mesh.boundary().size();
    reduce(havePatches, andOp<bool>());

    if (!havePatches)
    {
        // There are one or more processors without full complement of
        // patches

        DynamicList<polyPatch*> newPatches;

        if (haveMesh)   //Pstream::master())
        {
            forAll(mesh.boundary(), patchI)
            {
                newPatches.append
                (
                    mesh.boundaryMesh()[patchI].clone(mesh.boundaryMesh()).ptr()
                );
            }
        }
        else
        {
            forAll(patchEntries, patchI)
            {
                const entry& e = patchEntries[patchI];
                const word type(e.dict().lookup("type"));

                if
                (
                    type == processorPolyPatch::typeName
                 || type == processorCyclicPolyPatch::typeName
                )
                {}
                else
                {
                    dictionary patchDict(e.dict());
                    patchDict.set("nFaces", 0);
                    patchDict.set("startFace", 0);

                    newPatches.append
                    (
                        polyPatch::New
                        (
                            patchEntries[patchI].keyword(),
                            patchDict,
                            newPatches.size(),
                            mesh.boundaryMesh()
                        ).ptr()
                    );
                }
            }
        }
        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches);
    }

    //Pout<< "patches:" << endl;
    //forAll(mesh.boundary(), patchI)
    //{
    //    Pout<< "    type" << mesh.boundary()[patchI].type()
    //        << " size:" << mesh.boundary()[patchI].size()
    //        << " start:" << mesh.boundary()[patchI].start() << endl;
    //}
    //Pout<< endl;


    // Determine zones
    // ~~~~~~~~~~~~~~~

    wordList pointZoneNames(mesh.pointZones().names());
    Pstream::scatter(pointZoneNames);
    wordList faceZoneNames(mesh.faceZones().names());
    Pstream::scatter(faceZoneNames);
    wordList cellZoneNames(mesh.cellZones().names());
    Pstream::scatter(cellZoneNames);

    if (!haveMesh)
    {
        // Add the zones. Make sure to remove the old dummy ones first
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        List<pointZone*> pz(pointZoneNames.size());
        forAll(pointZoneNames, i)
        {
            pz[i] = new pointZone
            (
                pointZoneNames[i],
                labelList(0),
                i,
                mesh.pointZones()
            );
        }
        List<faceZone*> fz(faceZoneNames.size());
        forAll(faceZoneNames, i)
        {
            fz[i] = new faceZone
            (
                faceZoneNames[i],
                labelList(0),
                boolList(0),
                i,
                mesh.faceZones()
            );
        }
        List<cellZone*> cz(cellZoneNames.size());
        forAll(cellZoneNames, i)
        {
            cz[i] = new cellZone
            (
                cellZoneNames[i],
                labelList(0),
                i,
                mesh.cellZones()
            );
        }

        if (pz.size() && fz.size() && cz.size())
        {
            mesh.addZones(pz, fz, cz);
        }
    }

    return meshPtr;
}


// ************************************************************************* //

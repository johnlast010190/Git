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
    (c) 2019-2025 Engys Ltd.

Application
    helyxMap

Description
    Map two flow fields with arbitrary unstructured mesh.
    Must be run on maximum number of source and target processors.

\*---------------------------------------------------------------------------*/

#include "helyxMap.H"
#include "global/argList/argList.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "regionModel/regionProperties/regionProperties.H"
#include "redistributePar/loadOrCreateMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkProcessorFolders(const fileName& dir)
{
    const fileName lck(dir/"helyxMap.lock");

    while (Foam::isFile(lck))
    {
        Info<< "Found lock file " << lck << ". Waiting..." << endl;
        Foam::sleep(10);
    }

    label nProcs = 0;

    while (isDir(dir/"processor" + Foam::name(++nProcs)))
    {}

    // Wait for all procs to finish counting
    gMax(labelList(Pstream::nProcs(), 0));

    // Create dummy source folders
    if (nProcs < Pstream::nProcs())
    {
        instantList timeDirs;
        const word procFolder = "processor" + Foam::name(Pstream::myProcNo());

        if (Pstream::master())
        {
            const bool oldParRun = Pstream::parRun();

            // Only create lock file if it doesn't already exist
            if (!Foam::isFile(lck))
            {
                Info<< "Creating lock file" << endl;

                OFstream os(lck);
                os << "helyxMap: create dummy source folders\n";
                os.flush();
            }
            else
            {
                // Otherwise check again
                Foam::sleep(1);
                checkProcessorFolders(dir);
            }

            Pstream::parRun() = false;
            timeDirs = Time::findTimes(dir/procFolder, "constant");
            Pstream::parRun() = oldParRun;
        }

        Pstream::scatter(timeDirs);

        forAll(timeDirs, i)
        {
            mkDir(dir/procFolder/timeDirs[i].name());
        }
    }
    else if (nProcs > Pstream::nProcs())
    {
        FatalError
            << "Running with " << Pstream::nProcs() << " processors but case "
            << dir << " is decomposed in " << nProcs << " processors."
            << nl << "Please run the application with maximum number "
            << "of source and target processors."
            << exit(FatalError);
    }
}


autoPtr<fvMesh> detectAndGetMesh
(
    const word& regionName,
    const Time& runTime,
    boolList& haveMesh
)
{
    fileName meshSubDir;

    if (regionName == polyMesh::defaultRegion)
    {
        meshSubDir = polyMesh::meshSubDir;
    }
    else
    {
        meshSubDir = regionName/polyMesh::meshSubDir;
    }

    fileName masterInstDir;

    if (Pstream::master())
    {
        masterInstDir = runTime.findInstance
        (
            meshSubDir,
            "faces",
            IOobject::READ_IF_PRESENT
        );
    }
    Pstream::scatter(masterInstDir);

    // Check who has a mesh
    const fileName meshPath = runTime.path()/masterInstDir/meshSubDir/"faces";

    Info<< "    Checking for mesh in " << meshPath << nl << endl;

    haveMesh = boolList(Pstream::nProcs(), false);
    haveMesh[Pstream::myProcNo()] = isFile(meshPath);
    Pstream::allGatherList(haveMesh);

    Info<< "    Per processor mesh availability: " << haveMesh << endl;

    autoPtr<fvMesh> meshPtr =
        loadOrCreateMesh
        (
            IOobject
            (
                regionName,
                masterInstDir,
                runTime,
                Foam::IOobject::MUST_READ
            ),
            false    // removeDummyMesh
        );

    return meshPtr;
}


void writeDummyFields
(
    const boolList& haveMesh,
    const fvMesh& mesh,
    const Time& runTime
)
{
    autoPtr<fvMeshSubset> subsetterPtr;

    // Find last non-processor patch
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nonProcI = -1;

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI])) break;

        nonProcI++;
    }

    if (nonProcI == -1)
    {
        FatalErrorInFunction
            << "Cannot find non-processor patch on processor "
            << Pstream::myProcNo() << endl
            << " Current patches: " << patches.names()
            << abort(FatalError);
    }

    // Subset 0 cells, no parallel comms. This is used to create
    // zero-sized fields.
    subsetterPtr.reset(new fvMeshSubset(mesh));
    subsetterPtr().setLargeCellSubset(labelHashSet(0), nonProcI, false);

    PtrList<volScalarField> volScalarFields;
    PtrList<volVectorField> volVectorFields;

    IOobjectList objects(mesh, runTime.timeName());

    readFields(haveMesh, mesh, subsetterPtr, objects, volScalarFields);
    readFields(haveMesh, mesh, subsetterPtr, objects, volVectorFields);

    // Write empty field
    if (!haveMesh[Pstream::myProcNo()])
    {
        forAll(volScalarFields, i)
        {
            volScalarFields[i].write();
        }

        forAll(volVectorFields, i)
        {
            volVectorFields[i].write();
        }
    }

    subsetterPtr.reset();
}


void reportMemory(memInfo& mem, const word& desc, scalar& maxUsage)
{
    scalar memory = scalar(mem.update().size())*1.0e-3;

    if (memory > maxUsage)
    {
        scalar maxUsage = memory;

        Pout<< "Maximum memory usage now: " << maxUsage
            << " MB, " << desc << endl;
    }
}


int sourceMapSeq
(
    int argc,
    char *argv[],
    helyxMap& fieldMap,
    const Time& runTime,
    const fvMesh& mesh
)
{
    scalar maxUsage = 0;
    word desc = "before reading source mesh";
    memInfo mem;

    if (fieldMap.reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    if (fieldMap.reportMemUsage())
    {
        word desc = "after reading source mesh";
        reportMemory(mem, desc, maxUsage);
    }

    fieldMap.getFieldTypes(mesh, runTime.timeName());

    fieldMap.createSource(&mesh, &runTime);

    fieldMap.source().createFields(runTime.timeName());

    if (fieldMap.reportMemUsage())
    {
        word desc = "after creating source fields";
        reportMemory(mem, desc, maxUsage);
    }

    fieldMap.source().setInputs();

    fieldMap.getWallDistPatchs(mesh);

    const volScalarField y =
        fieldMap.wallDistPatchs().size() == 0
      ? wallDist::New(mesh).y()
      : wallDist::New(mesh, fieldMap.wallDistPatchs(), "wall").y();

    if (fieldMap.reportMemUsage())
    {
        word desc = "after getting wall distance in source";
        reportMemory(mem, desc, maxUsage);
    }

    fieldMap.source().buildKdTrees(y);

    if (fieldMap.reportMemUsage())
    {
        word desc = "after building kdTree";
        reportMemory(mem, desc, maxUsage);
    }

    fieldMap.source().storeFields();

    if (fieldMap.reportMemUsage())
    {
        word desc = "after store fields";
        reportMemory(mem, desc, maxUsage);
    }

    return 0;
}


int sourceMap
(
    int argc,
    char *argv[],
    helyxMap& fieldMap,
    const Time& runTime,
    const fvMesh& mesh,
    boolList& haveMesh
)
{
    Info<< "Reading source" << endl;

    scalar maxUsage = 0;
    word desc = "before reading source mesh";
    memInfo mem;

    if (fieldMap.reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    const bool allHaveMesh = (findIndex(haveMesh, false) == -1);
    if (!allHaveMesh)
    {
        writeDummyFields(haveMesh, mesh, runTime);
    }

    fieldMap.getFieldTypes(mesh, runTime.timeName());

    fieldMap.createSource(&mesh, &runTime);

    fieldMap.source().createFields(runTime.timeName());

    fieldMap.source().setInputs();

    if (fieldMap.reportMemUsage())
    {
        desc = "after reading source fields";
        reportMemory(mem, desc, maxUsage);
    }

    if (fieldMap.wDistMap_)
    {
        fieldMap.getWallDistPatchs(mesh);

        fieldMap.wait();

        const volScalarField y =
            fieldMap.wallDistPatchs().size() == 0
          ? wallDist::New(mesh).y()
          : wallDist::New(mesh, fieldMap.wallDistPatchs(), "wall").y();

        scalar maxWDist = gMax(y);

        fieldMap.source().setWallDist(y, maxWDist);
    }

    if (fieldMap.reportMemUsage())
    {
        desc = "after calculating wall distance in source case";
        reportMemory(mem, desc, maxUsage);
    }

    fieldMap.source().storeFields();

    if (fieldMap.reportMemUsage())
    {
        desc = "before exiting source case";
        reportMemory(mem, desc, maxUsage);
    }

    if (!haveMesh[Pstream::myProcNo()])
    {
        // We created a dummy mesh file above. Delete it.
        const fileName meshFiles = runTime.path();
        rmDir(meshFiles);
    }

    if (!allHaveMesh)
    {
        // Wait for all procs to finish
        gMax(labelList(Pstream::nProcs(), 0));

        if (Pstream::master())
        {
            rm(runTime.rootPath()/"helyxMap.lock");
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "sourceCase",
        "dir",
        "Specify source case directory"
    );
    argList::addOption
    (
        "mapFixedBCs",
        "List<word>",
        "Specify the boundary conditions to map from"
    );
    argList::addOption
    (
        "sourceTime",
        "scalar|latestTime",
        "Specify the time step/iteration of the source case to map from \
            (default: latestTime)"
    );
    argList::addOption
    (
        "targetTime",
        "scalar",
        "Specify the time step/iteration of the target case to map to \
            (default: 0)"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "Specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "Specify the target region"
    );
    argList::addBoolOption
    (
        "allRegions",
        "Operate on all regions in target regionProperties"
    );
    argList::addOption
    (
        "regionMaps",
        "List<Pair<word>>",
        "Specify source region to target region maps"
    );
    argList::addOption
    (
        "mapScalarFields",
        "list",
        "Specify the scalar field names to be mapped, e.g. '(p nut)'"
    );
    argList::addOption
    (
        "mapVectorFields",
        "list",
        "Specify the vector field names to be mapped, e.g. '(U)'"
    );
    argList::addOption
    (
        "nwdist",
        "label",
        "Specify the number of divisions in the wall distance range[0,1]"
    );
    argList::addOption
    (
        "rhoRefSource",
        "scalar",
        "Specify the reference density value of the source field"
    );
    argList::addOption
    (
        "rhoRefTarget",
        "scalar",
        "Specify the reference density value of the target field"
    );
    argList::addOption
    (
        "UrefSource",
        "vector",
        "Specify the reference velocity of the source field, e.g. '(1 0 0)'"
    );
    argList::addOption
    (
        "UrefTarget",
        "vector",
        "Specify the reference velocity of the target field, e.g. '(1 0 0)'"
    );
    argList::addOption
    (
        "UrotDegreeFromSource",
        "scalar",
        "Specify the in-flow velocity vector rotate degree \
            from source to target"
    );
    argList::addOption
    (
        "fieldTypes",
        "HashTable<word>",
        "Specify the property types of each maping field, \
            e.g. '(p pressure U velocity k turbEnergy)'"
    );
    argList::addOption
    (
        "alphaMax",
        "scalar",
        "Specify the upper wall-distance value to be included in the \
            wall-distance-based search trees"
    );
    argList::addBoolOption
    (
        "mapBoundary",
        "Specify whether the map includes boundaries"
    );
    argList::addBoolOption
    (
        "interpolation",
        "Specify whether interpolation will be used in the mapping"
    );
    argList::addOption
    (
        "function",
        "map|validate",
        "Specify the purpose of running helyxMap:\n\
            'map' is for mapping two flows or map a flow onto a grid,\n\
            'validate' validates the method by mapping a field into itself."
    );
    argList::addNote
    (
        "helyxMap maps an arbitrary number of flow fields from an unstructed\n\
            mesh domain onto another unstructured domain. The two domains\n\
            must be similar in nature but do not have to be of the same size.\n\
            The boundary conditions do not need to be the same. The in-flow\n\
            velocity can be rotated and scaled from the source domain\n\
            to the target domain."
    );

    argList::noCheckProcessorDirectories();

    #if !defined( WIN32 ) && !defined( WIN64 )
    #include "include/addProfilingOption.H"
    #endif

    #include "include/setRootCase.H"

    // Check target processors
    checkProcessorFolders(args.rootPath()/args.globalCaseName());

    Foam::Time runTimeTarget(Foam::Time::controlDictName, args);

    helyxMap fieldMap;

    IOdictionary dict
    (
        IOobject
        (
            "helyxMapDict",
            runTimeTarget.system(),
            "",
            runTimeTarget,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fieldMap.setInput(dict, args);

    fieldMap.myId_     = Pstream::myProcNo();
    fieldMap.masterId_ = Pstream::masterNo();
    fieldMap.nProcs_   = Pstream::nProcs();

    if (Pstream::parRun())
    {
        fieldMap.parInit(fieldMap.nProcs_);
    }

    fieldMap.setOptions(args);

    fieldMap.mapCase_ = cwd();

    const fileName rootDirSource = fileName(fieldMap.sourceCase_).toAbsolute();

    if (!isDir(rootDirSource.path()))
    {
        FatalErrorInFunction
            << "Source case directory: " << rootDirSource << " does not exist!"
            << exit(FatalError);
    }

    fileName caseName = "./";
    if (Pstream::parRun())
    {
        caseName = caseName/fileName(word("processor") + name(fieldMap.myId_));
    }

    Info<< "sourceCase: " << rootDirSource << endl;

    // Check sourceCase processors
    checkProcessorFolders(rootDirSource);

    Time runTimeSource(Time::controlDictName, rootDirSource, caseName);

    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();

    if (fieldMap.mapTimeName_ == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(fieldMap.mapTimeName_);

        sourceTimeIndex =
            Time::findClosestTimeIndex(sourceTimes, readScalar(is));
    }

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    fieldMap.mapTimeName_ = runTimeSource.timeName();

    Info<< "sourceTime: " << runTimeSource.timeName() << endl;

    const word oldTargetTime = runTimeTarget.timeName();

    instantList targetTimes = runTimeTarget.times();
    label targetTimeIndex = runTimeTarget.timeIndex();

    if (fieldMap.tgtTimeName_ == "latestTime")
    {
        targetTimeIndex = targetTimes.size() - 1;

        runTimeTarget.setTime
        (
            targetTimes[targetTimeIndex],
            targetTimes[targetTimeIndex].value()
        );
    }
    else
    {
        IStringStream is(fieldMap.tgtTimeName_);

        const instant targetTime(readScalar(is), fieldMap.tgtTimeName_);

        runTimeTarget.setTime(targetTime, 0);
    }

    fieldMap.tgtTimeName_ = runTimeTarget.timeName();

    if
    (
        oldTargetTime != fieldMap.tgtTimeName_
    && !isDir(fileName(args.path()/fieldMap.tgtTimeName_))
    )
    {
        if
        (
           !cp
            (
                fileName(args.path()/oldTargetTime),
                fileName(args.path()/fieldMap.tgtTimeName_),
                false    // followLink
            )
        )
        {
            FatalErrorInFunction
                << "Cannot copy time directory " << oldTargetTime
                << " to " << fieldMap.tgtTimeName_
                << exit(FatalError);
        }
    }

    Info<< "targetTime: " << runTimeTarget.timeName() << endl;

    // Warn the user if targetTime is different from startTime in controlDict
    if
    (
        Foam::name(runTimeTarget.startTime().value())
     != runTimeTarget.timeName()
    )
    {
        Warning
            << "Mapping to time " << runTimeTarget.timeName()
            << " which is different from startTime "
            << runTimeTarget.startTime().value() << endl;
    }

    if (fieldMap.allRegions_)
    {
        Info<< "\nMapping for all regions in target regionProperties" << endl;

        fieldMap.sourceRegions_.resize(0);
        fieldMap.targetRegions_.resize(0);

        regionProperties rp(runTimeTarget);

        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, regioni)
            {
                if (findIndex(fieldMap.targetRegions_, regions[regioni]) == -1)
                {
                    fieldMap.targetRegions_.append(regions[regioni]);
                    fieldMap.sourceRegions_.append(regions[regioni]);
                }
            }
        }

        fieldMap.targetRegionDirs_ = fieldMap.targetRegions_;
    }

    scalar maxUsage = 0;
    word desc = "before reading target mesh";
    memInfo mem;

    forAll(fieldMap.sourceRegions_, regioni)
    {
        const word& sourceRegion = fieldMap.sourceRegions_[regioni];
        const word& targetRegion = fieldMap.targetRegions_[regioni];

        Info<< "\nMapping source region: "<< sourceRegion
            << "\n     to target region: "<< targetRegion << endl;

        autoPtr<fvMesh> sourceMeshPtr;

        if (Pstream::parRun())
        {
            boolList haveMesh;

            sourceMeshPtr =
                detectAndGetMesh(sourceRegion, runTimeSource, haveMesh);

            sourceMap
            (
                argc,
                argv,
                fieldMap,
                runTimeSource,
                sourceMeshPtr(),
                haveMesh
            );

            fieldMap.wait();
        }
        else
        {
            chDir(fieldMap.sourceCase_);

            sourceMeshPtr.reset
            (
                new fvMesh
                (
                    IOobject
                    (
                        sourceRegion,
                        runTimeSource.timeName(),
                        runTimeSource,
                        IOobject::MUST_READ
                    )
                )
            );

            sourceMapSeq(argc, argv, fieldMap, runTimeSource, sourceMeshPtr());
        }

        Info<< "\nReading target" << endl;

        boolList targetHaveMesh;

        autoPtr<fvMesh> targetMeshPtr =
            detectAndGetMesh(targetRegion, runTimeTarget, targetHaveMesh);

        fvMesh& targetMesh = targetMeshPtr();

        if (fieldMap.reportMemUsage())
        {
            desc = "after reading target mesh";
            reportMemory(mem, desc, maxUsage);
        }

        const bool allHaveMesh = (findIndex(targetHaveMesh, false) == -1);
        if (!allHaveMesh)
        {
            writeDummyFields(targetHaveMesh, targetMesh, runTimeTarget);
        }

        fieldMap.createTarget(&targetMesh, &runTimeTarget);

        fieldMap.target().createFields(fieldMap.tgtTimeName_);
        Info<< "\nTarget field created." << nl << endl;

        fieldMap.checkFieldDimensions();

        sourceMeshPtr.reset();

        if (fieldMap.reportMemUsage())
        {
            desc = "after reading target fields";
            reportMemory(mem, desc, maxUsage);
        }

        fieldMap.target().setInputs();

        if (fieldMap.scaleSource())
        {
            fieldMap.source().transformPoints
            (
                fieldMap.target().gBox(),
                fieldMap.source().xyz()
            );

            Info<< "Points transformed..." << endl;
        }

        // Build searchTree for each processor
        fieldMap.getDomainFields();

        if (fieldMap.reportMemUsage())
        {
            desc = "after get domain fields";
            reportMemory(mem, desc, maxUsage);
        }

        Info<< "\nDomain field obtained..." << nl << endl;

        if (fieldMap.wDistMap_)
        {
            fieldMap.getWallDistPatchs(targetMesh);

            if (fieldMap.wallDistPatchs().size() == 0)
            {
                Info<< "Setting wall distance..." << endl;

                const volScalarField y = wallDist::New(targetMesh).y();
                scalar maxWDist = gMax(y);
                fieldMap.target().setWallDist(y, maxWDist);

                if (Pstream::parRun())
                {
                    fieldMap.wait();
                }
            }
            else
            {
                const volScalarField y =
                    wallDist::New
                    (
                        targetMesh,
                        fieldMap.wallDistPatchs(),
                        "wall"
                    ).y();

                scalar maxWDist = gMax(y);
                fieldMap.target().setWallDist(y, maxWDist);

                if (Pstream::parRun())
                {
                    fieldMap.wait();
                }
            }
        }

        if (fieldMap.reportMemUsage())
        {
            desc = "after creating wall distance field in target";
            reportMemory(mem, desc, maxUsage);
        }

        if (Pstream::parRun())
        {
            Info<< "Building parallel search trees" << endl;
            fieldMap.buildKdTreesPar();

            if (fieldMap.reportMemUsage())
            {
                desc = "after building parallel search trees in target";
                reportMemory(mem, desc, maxUsage);
            }

            // Parallel-run map fields
            fieldMap.mapFields();

            if (fieldMap.reportMemUsage())
            {
                desc = "after mapping fields in target";
                reportMemory(mem, desc, maxUsage);
            }
        }
        else
        {
            if (fieldMap.function_ == "validate")
            {
                Info<< "Validating helyxMap" << endl;

                fieldMap.searchError();
                fieldMap.mapError("p");
                fieldMap.mapError("nuTilda");
                fieldMap.UmapError("U");
            }
            else
            {
                Info<< "Mapping fields" << nl << endl;
                fieldMap.mapFields();

                if (fieldMap.reportMemUsage())
                {
                    desc = "after mapping fields in target";
                    reportMemory(mem, desc, maxUsage);
                }
            }
        }

        runTimeTarget.writeNow();

        fieldMap.clearSourceTarget();

        if (!targetHaveMesh[Pstream::myProcNo()])
        {
            // We created a dummy mesh file above. Delete it.
            const fileName meshFiles = runTimeTarget.path();
            rmDir(meshFiles);
        }

        targetMeshPtr.reset();
    }

    Info<< "ExecutionTime = " << runTimeTarget.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTimeTarget.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

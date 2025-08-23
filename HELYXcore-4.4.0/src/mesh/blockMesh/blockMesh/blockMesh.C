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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "blockMesh/blockMesh.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicTransform.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "coordinate/systems/coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitch(blockMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMesh::blockMesh
(
   const IOdictionary& dict,
   const word& regionName,
   const bool readBoundaryFile
)
:
    meshDict_(dict),
    verboseOutput(meshDict_.lookupOrDefault<Switch>("verbose", true)),
    geometry_
    (
        IOobject
        (
            "geometry",                 // dummy name
            meshDict_.time().constant(),     // instance
            "geometry",                 // local
            meshDict_.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshDict_.found("geometry")
      ? meshDict_.subDict("geometry")
      : dictionary(),
        true
    ),
    scaleFactor_(1.0),
    blockVertices_
    (
        meshDict_.lookup("vertices"),
        blockVertex::iNew(meshDict_, geometry_)
    ),
    vertices_(Foam::vertices(blockVertices_)),
    topologyPtr_(createTopology(meshDict_, regionName, readBoundaryFile))
{
    Switch fastMerge(meshDict_.lookupOrDefault<Switch>("fastMerge", false));

    if (fastMerge)
    {
        calcMergeInfoFast();
    }
    else
    {
        calcMergeInfo();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMesh::~blockMesh()
{
    delete topologyPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMesh::verbose(const bool on)
{
    verboseOutput = on;
}


const Foam::pointField& Foam::blockMesh::vertices() const
{
    return vertices_;
}


const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorInFunction
            << "topologyPtr_ not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::PtrList<Foam::dictionary> Foam::blockMesh::patchDicts() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();

    PtrList<dictionary> patchDicts(patchTopologies.size());

    forAll(patchTopologies, patchi)
    {
        autoPtr<polyPatch> ppPtr =
            patchTopologies[patchi].clone(topology().boundaryMesh());

        if (isA<cyclicTransform>(ppPtr()))
        {
            refCast<cyclicTransform>(ppPtr()) =
                transformer::scaling(scaleFactor_*tensor::I)
              & refCast<cyclicTransform>(ppPtr());
        }

        OStringStream os;
        ppPtr->write(os);
        IStringStream is(os.str());
        patchDicts.set(patchi, new dictionary(is));
    }

    return patchDicts;
}


Foam::scalar Foam::blockMesh::scaleFactor() const
{
    return scaleFactor_;
}


const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();

        // Use reference frames
        if (meshDict_.found("referenceFrame"))
        {
            const word frameName = meshDict_.lookup<word>("referenceFrame");
            Info<<"Using reference frame " << frameName <<endl;

            autoPtr<dictionary> meshObjectDict;
            const word dictName = "meshObjects";
            autoPtr<IOobject> meshObjectIO
            (
                new IOobject
                (
                    dictName,
                    meshDict_.db().time().caseSystem(),
                    meshDict_.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            if (!meshObjectIO->typeHeaderOk<localIOdictionary>())
            {
                // Fallback: look for definition at the top level
                autoPtr<IOobject> defaultRegionIO
                (
                    new IOobject
                    (
                        dictName,
                        meshDict_.db().time().caseSystem(),
                        meshDict_.db().time(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     )
                 );
                if (defaultRegionIO->typeHeaderOk<localIOdictionary>())
                {
                    meshObjectIO = defaultRegionIO;
                }
            }

            // Use objectPath rather than filePath so a useful error is given
            // if file is not found
            IFstream is(meshObjectIO->objectPath());

            meshObjectDict.reset(new dictionary(is));

            const coordinateSystem frame = frameData::frameToGlobal
            (
                meshObjectDict().subDict(frameName),
                meshDict_.time()
            );
            points_ = frame.transformPoint(points_);
        }
    }

    return points_;
}

const Foam::cellShapeList& Foam::blockMesh::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::faceListList& Foam::blockMesh::patches() const
{
    if (patches_.empty())
    {
        createPatches();
    }

    return patches_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    return topology().boundaryMesh().names();
}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    label num = 0;

    forAll(*this, blocki)
    {
        if (operator[](blocki).zoneName().size())
        {
            num++;
        }
    }

    return num;
}


void Foam::blockMesh::writeTopology(Ostream& os) const
{
    const pointField& pts = topology().points();

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    const edgeList& edges = topology().edges();

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        os << "l " << e.start() + 1 << ' ' << e.end() + 1 << endl;
    }
}

// ************************************************************************* //

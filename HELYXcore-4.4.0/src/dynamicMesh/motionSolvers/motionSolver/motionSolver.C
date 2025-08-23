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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/motionSolverList/motionSolverList.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "twoDPointCorrector/twoDPointCorrector.H"
#include "meshes/polyMesh/zones/ZoneMesh/cellZoneMesh.H"
#include "sets/topoSets/cellSet.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    mesh_(mesh),
    coeffDict_(dict.optionalSubDict(type + "Coeffs")),
    moveAllPoints_(false)
{
    updateSetPointIndices();
    initTransforms();
}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::clone() const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh,
    const dictionary& solverDict
)
{
    const word solverName
    (
        solverDict.found("motionSolver")
      ? solverDict.lookup("motionSolver")
      : solverDict.lookup("solver")
    );

    Info<< "Selecting motion solver: " << solverName << endl;

    libs.open(solverDict, "motionSolverLibs");

    const auto ctor =
        ctorTableLookup
        (
            "solver type",
            dictionaryConstructorTable_(),
            solverName
        );

    return autoPtr<motionSolver>(ctor(mesh, solverDict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolver::~motionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::motionSolver::initTransforms()
{
    isSolidBody_ = false;

    if
    (
        coeffDict().found("referenceFrame")
     || coeffDict().found(solidBodyMotionFunction::typeName)
    )
    {
        if (coeffDict().found(solidBodyMotionFunction::typeName))
        {
            SBMFs_.setSize(1);
            SBMFs_.set
            (
                0,
                solidBodyMotionFunction::New(coeffDict(), mesh().time())
            );
            isSolidBody_ = true;

            DeprecationWarningInFunction
            (
                solidBodyMotionFunction::typeName,
                "mesh motion type",
                40200,
                "Please replace it by using referenceFrame."
            );
        }
        else
        {
            refFrames_.setSize(1);
            refFrames_.set
            (
                0,
                coordinateFrame::lookupNew
                (
                    dynamic_cast<const fvMesh&>(mesh()),
                    coeffDict()
                )
            );

            refFrames_[0].resetDynamic(true);

            // Future check for dynamic motion
            if (!refFrames_[0].anyDynamic())
            {
                FatalErrorInFunction
                    << "coordinateFrame " << refFrames_[0].name()
                    << " is linked to dynamic mesh but it is not dynamic."
                    << abort(FatalError);
            }
        }
    }
    else
    {
        wordList names;

        if (coeffDict().found("referenceFrames"))
        {
            names = coeffDict().lookup<wordList>("referenceFrames");
            refFrames_.setSize(names.size());
        }
        else
        {
            forAllConstIter(dictionary, coeffDict(), iter)
            {
                const word& name = iter().keyword();
                if (iter().isDict() && name != word("solver"))
                {
                    if
                    (
                        coeffDict().subDict(name).found
                        (
                            solidBodyMotionFunction::typeName
                        )
                    )
                    {
                        names.append(name);
                    }
                }
            }
            SBMFs_.setSize(names.size());
            isSolidBody_ = true;
            if (names.size())
            {
                DeprecationWarningInFunction
                (
                    solidBodyMotionFunction::typeName,
                    "mesh motion type",
                    40200,
                    "Please replace it by using referenceFrame."
                );
            }
        }
        forAll(names, namei)
        {
            const word& name = names[namei];

            if (isSolidBody_)
            {
                SBMFs_.set
                (
                    namei,
                    solidBodyMotionFunction::New
                    (
                        coeffDict().subDict(name),
                        mesh().time()
                    )
                );
            }
            else
            {
                refFrames_.set
                (
                    namei,
                    &coordinateFrame::New
                    (
                        dynamic_cast<const fvMesh&>(mesh()),
                        name
                    )
                );
                refFrames_[namei].resetDynamic(true);

                // Future check for dynamic motion
                if (!refFrames_[namei].anyDynamic())
                {
                    FatalErrorInFunction
                        << "coordinateFrame " << refFrames_[namei].name()
                        << " is linked to dynamic mesh but it is not dynamic."
                        << abort(FatalError);
                }
            }
        }
    }
}


bool Foam::motionSolver::isIncrementalMotion() const
{
    // Please note all motions have to be relative or absolute
    if (SBMFs_.size())
    {
        return SBMFs_[0].isIncrementalMotion();
    }
    else if (refFrames_.size())
    {
        return refFrames_[0].isIncrementalMotion();
    }
    return false;
}


Foam::septernion Foam::motionSolver::transformation
(
    const label transformationI,
    const wordList& nestedFrames
) const
{
    if (isSolidBody_)
    {
        return SBMFs_[transformationI].transformation();
    }

    // Necessary to support old solvers
    // (helyxSolve has it's own updateStates)
    coordinateFrame::updateStates(mesh().thisDb());

    return refFrames_[transformationI].transformation();
}


Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();
    return curPoints();
}


void Foam::motionSolver::twoDCorrectPoints(pointField& p) const
{
    twoDPointCorrector::New(mesh_).correctPoints(p);
}


void Foam::motionSolver::updateSetPointIndices()
{
    const word cellZoneName =
        coeffDict().lookupOrDefault<word>("cellZone", "none");

    const word cellSetName =
        coeffDict().lookupOrDefault<word>("cellSet", "none");

    if ((cellZoneName != "none") && (cellSetName != "none"))
    {
        FatalIOErrorInFunction(coeffDict())
            << "Either cellZone OR cellSet can be supplied, but not both. "
            << "If neither is supplied, all cells will be included"
            << exit(FatalIOError);
    }

    labelList cellIDs;
    if (cellZoneName != "none")
    {
        const label zoneID = mesh_.cellZones().findZoneID(cellZoneName);

        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Unable to find cellZone " << cellZoneName
                << ".  Valid cellZones are:"
                << mesh_.cellZones().names()
                << exit(FatalError);
        }

        cellIDs = mesh_.cellZones()[zoneID];
    }

    if (cellSetName != "none")
    {
        cellSet set(mesh_, cellSetName);
        cellIDs = set.toc();
    }

    const label nCells = returnReduce(cellIDs.size(), sumOp<label>());
    moveAllPoints_ = nCells == 0;

    if (moveAllPoints_)
    {
        setPointIndices_ = identity(mesh_.nPoints());
    }
    else
    {
        boolList movePts(mesh_.nPoints(), false);

        forAll(cellIDs, i)
        {
            label celli = cellIDs[i];
            const cell& c = mesh_.cells()[celli];
            forAll(c, j)
            {
                const face& f = mesh_.faces()[c[j]];
                forAll(f, k)
                {
                    movePts[f[k]] = true;
                }
            }
        }

        syncTools::syncPointList(mesh_, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(mesh_.nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        setPointIndices_.transfer(ptIDs);
    }
}


void Foam::motionSolver::topoChange(const polyTopoChangeMap& map)
{}


// ************************************************************************* //

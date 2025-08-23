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
    (c) 2019 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/motionSolverList/motionSolverList.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolverList, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        motionSolverList,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolverList::motionSolverList
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName)
{
    const dictionary& solversDict = dict.subDict("solvers");

    motionSolvers_.setSize(solversDict.size());

    label zonei = 0;
    forAllConstIter(dictionary, solversDict, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            motionSolvers_.set(zonei, motionSolver::New(mesh, subDict).ptr());

            zonei++;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolverList::~motionSolverList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolverList::curPoints() const
{
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(mesh().nPoints(), Zero);

        forAll(motionSolvers_, i)
        {
            disp += motionSolvers_[i].curPoints() - mesh().points();
        }

        return mesh().points() + disp;
    }
    else
    {
        return mesh().points();
    }
}


void Foam::motionSolverList::solve()
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].solve();
    }
}


void Foam::motionSolverList::topoChange(const polyTopoChangeMap& map)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].topoChange(map);
    }
}


void Foam::motionSolverList::mapMesh(const polyMeshMap& map)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].mapMesh(map);
    }
}


void Foam::motionSolverList::distribute
(
    const polyDistributionMap& map
)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].distribute(map);
    }
}


void Foam::motionSolverList::movePoints(const pointField& points)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].movePoints(points);
    }
}


bool Foam::motionSolverList::write() const
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].write();
    }

    return motionSolver::write();
}


// ************************************************************************* //

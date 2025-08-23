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
    (c) 2014 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.H"

namespace Foam
{

defineTypeNameAndDebug(meshOptimizationMethod, 0);

defineRunTimeSelectionTable(meshOptimizationMethod, dictionary);

/*---------------------------------------------------------------------------*\
                   Class meshOptimizationMethod Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
labelList meshOptimizationMethod::compactList(const boolList& bL)
{
    label num = 0;
    forAll(bL, i)
    {
        if (bL[i])
        {
            num++;
        }
    }

    labelList set(num);
    label counter = 0;
    forAll(bL, i)
    {
        if (bL[i])
        {
            set[counter] = i;
            counter++;
        }
    }
    return set;
}

void meshOptimizationMethod::incrementActivePointsSet(boolList& activePoints)
{
    boolList activeCells(mesh_.cells().size(), false);

    forAll(activePoints, pI)
    {
        if (activePoints[pI])
        {
            forAll(mesh_.pointCells()[pI], cI)
            {
                const label& cL = mesh_.pointCells()[pI][cI];
                activeCells[cL] = true;
            }
        }
    }

    forAll(mesh_.cells(), cI)
    {
        if (activeCells[cI])
        {
            forAll(mesh_.cellPoints()[cI], pI)
            {
                const label& pL = mesh_.cellPoints()[cI][pI];
                activePoints[pL] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        activePoints,
        maxEqOp<bool>(),
        false
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshOptimizationMethod::meshOptimizationMethod
(
    const dictionary& dict,
    const fvMesh& mesh,
    const meshMetric& metric
)
:
    mesh_(mesh)
{}

Foam::meshOptimizationMethod::~meshOptimizationMethod()
{}

const activeSet& Foam::meshOptimizationMethod::getActiveSet() const
{
    return activeSet_;
}

} /* namespace Foam */

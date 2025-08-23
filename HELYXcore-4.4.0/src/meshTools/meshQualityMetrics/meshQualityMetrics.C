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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshQualityMetrics/meshQualityMetrics.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshTools/meshTools.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshQualityMetrics, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshQualityMetrics::readMetricThresholds()
{
    minVol_ = mqDict_.lookupOrDefault<scalar>("minVol", 1e-16);
    minArea_ = mqDict_.lookupOrDefault<scalar>("minArea", 2e-13);
    maxNonOrtho_ = mqDict_.lookupOrDefault<scalar>("maxNonOrtho", 70);
    maxSkew_ = mqDict_.lookupOrDefault<scalar>("maxSkew", 100);
    minDeterminant_ = mqDict_.lookupOrDefault<scalar>("minDeterminant", 1e-4);
    minFaceWeight_ = mqDict_.lookupOrDefault<scalar>("minFaceWeight", 0.08);
    reportStatistics_ = mqDict_.lookupOrDefault<Switch>
    (
        "reportStatistics", false
    );
    writeField_ = mqDict_.lookupOrDefault<Switch>
    (
        "writeField", false
    );
}


void Foam::meshQualityMetrics::calcMetrics() const
{
    badQualityCellsPtr_ = new boolList(mesh_.nCells(), false);
    boolList& badQualityCells = *badQualityCellsPtr_;

    cellSet errorCells(mesh_, "errorCells", mesh_.nCells()/100+1);
    faceSet errorFaces(mesh_, "errorFaces", mesh_.nFaces()/100+1);

    if (reportStatistics_) Info<< this->name() << tab << ":" << endl;

    //- cell volume
    mesh_.checkCellVolumes(reportStatistics_, &errorCells);

    //- face areas
    mesh_.checkFaceAreas(reportStatistics_, minArea_, &errorFaces);

    //- face orthogonality
    mesh_.checkFaceOrthogonality(reportStatistics_, &errorFaces, maxNonOrtho_);

    //- face skewness
    mesh_.checkFaceSkewness(reportStatistics_, &errorFaces, maxSkew_);

    //- face pyramids
    mesh_.checkFacePyramids(reportStatistics_, minVol_, &errorFaces);

    //- Cell determinant
    mesh_.checkCellDeterminant(reportStatistics_, &errorFaces, minDeterminant_);

    //- Face weights
    mesh_.checkFaceWeight(reportStatistics_, minFaceWeight_, &errorFaces);

    //- mark the own-nei of poor quality faces
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    forAllConstIter(labelHashSet, errorFaces, iter)
    {
        label fII = iter.key();
        badQualityCells[own[fII]] = true;
        if (fII<mesh_.nInternalFaces())
        {
            badQualityCells[nei[fII]] = true;
        }
    }

    forAll(errorCells, cI)
    {
        badQualityCells[errorCells[cI]] = true;
    }
    if (reportStatistics_)
    {
        label nBadCells(0);
        label nCells(this->mesh().nCells());
        forAll(badQualityCells, cI)
        {
            if (badQualityCells[cI]) nBadCells++;
        }
        reduce(nBadCells, sumOp<label>());
        reduce(nCells, sumOp<label>());
        if (nCells>0)
        {
            Info<< tab << "Marked cells with poor quality "
                 << nBadCells << " out of " << nCells << " : "
                 << ((scalar(nBadCells)/scalar(nCells))*100) << "%" << endl;
        }
    }
}


void Foam::meshQualityMetrics::clearMetrics() const
{
    deleteDemandDrivenData(badQualityCellsPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshQualityMetrics::meshQualityMetrics
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& mqName
)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, meshQualityMetrics>(mesh, mqName),
    badQualityCellsPtr_(nullptr),
    mqDict_(dict),
    reportStatistics_(false),
    writeField_(false),
    minVol_(Zero),
    minArea_(Zero),
    maxNonOrtho_(Zero),
    maxSkew_(Zero),
    minDeterminant_(Zero),
    minFaceWeight_(Zero)
{
    readMetricThresholds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshQualityMetrics::~meshQualityMetrics()
{
    clearMetrics();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::meshQualityMetrics::topoChange(const polyTopoChangeMap&)
{
    calcMetrics();
}


void Foam::meshQualityMetrics::mapMesh(const polyMeshMap&)
{
    calcMetrics();
}


void Foam::meshQualityMetrics::distribute(const polyDistributionMap&)
{
    calcMetrics();
}


bool Foam::meshQualityMetrics::movePoints()
{
    return true;
}


const Foam::boolList& Foam::meshQualityMetrics::badQualityCells() const
{
    if (!badQualityCellsPtr_)
    {
        calcMetrics();
    }

    return *badQualityCellsPtr_;
}

// ************************************************************************* //

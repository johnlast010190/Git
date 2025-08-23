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
    (c) 2010-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/pointVolInterpolation/pointVolInterpolation.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "include/demandDrivenData.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "fields/Fields/Field/SubField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointVolInterpolation, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointVolInterpolation::makeWeights() const
{
    if (volWeightsPtr_)
    {
        FatalErrorInFunction
            << "weighting factors already calculated"
            << abort(FatalError);
    }
    if (debug)
    {
        Info<< "pointVolInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const pointField& points = vMesh().points();
    const labelListList& cellPoints = vMesh().cellPoints();
    const vectorField& cellCentres = vMesh().cellCentres();

    // Allocate storage for weighting factors
    volWeightsPtr_ = new FieldField<Field, scalar>(cellCentres.size());
    FieldField<Field, scalar>& weightingFactors = *volWeightsPtr_;

    forAll(weightingFactors, pointi)
    {
        weightingFactors.set
        (
            pointi,
            new scalarField(cellPoints[pointi].size())
        );
    }

    // Calculate inverse distances between points and cell centres
    // and store in weighting factor array
    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        forAll(curCellPoints, cellPointI)
        {
            weightingFactors[cellI][cellPointI] = 1.0/
                mag
                (
                    cellCentres[cellI] - points[curCellPoints[cellPointI]]
                );
        }
    }
    scalarField pointVolSumWeights(cellCentres.size(), Zero);

    // Calculate weighting factors using inverse distance weighting
    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        forAll(curCellPoints, cellPointI)
        {
            pointVolSumWeights[cellI] += weightingFactors[cellI][cellPointI];
        }
    }
    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        forAll(curCellPoints, cellPointI)
        {
            weightingFactors[cellI][cellPointI] /= pointVolSumWeights[cellI];
        }
    }

    if (debug)
    {
        Info<< "pointVolInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// Do what is necessary if the mesh has changed topology
void Foam::pointVolInterpolation::clearAddressing() const
{
    deleteDemandDrivenData(patchInterpolatorsPtr_);
}


// Do what is necessary if the mesh has moved
void Foam::pointVolInterpolation::clearGeom() const
{
    deleteDemandDrivenData(volWeightsPtr_);
}


const Foam::PtrList<Foam::primitivePatchInterpolation>&
Foam::pointVolInterpolation::patchInterpolators() const
{
    if (!patchInterpolatorsPtr_)
    {
        const fvBoundaryMesh& bdry = vMesh().boundary();
        patchInterpolatorsPtr_ =
            new PtrList<primitivePatchInterpolation>(bdry.size());

        forAll(bdry, patchI)
        {
            const polyPatch& pp = bdry[patchI].patch();
            if (isA<directPolyPatch>(pp))
            {
                const directPolyPatch& dpp =
                    refCast<const directPolyPatch>(pp);

                patchInterpolatorsPtr_->set
                (
                    patchI,
                    new primitivePatchInterpolation(dpp)
                );
            }
        }
    }
    return *patchInterpolatorsPtr_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::pointVolInterpolation::pointVolInterpolation
(
    const fvMesh& vm
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, pointVolInterpolation>(vm),
    pointMesh_(pointMesh::New(vm)),
    fvMesh_(vm),
    volWeightsPtr_(nullptr),
    patchInterpolatorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::pointVolInterpolation::~pointVolInterpolation()
{
    clearAddressing();
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return point weights
const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::pointVolInterpolation::volWeights() const
{
    // If weighting factor array not present then construct
    if (!volWeightsPtr_)
    {
        makeWeights();
    }

    return *volWeightsPtr_;
}

// Do what is necessary if the mesh has moved
bool Foam::pointVolInterpolation::movePoints()
{
    clearGeom();
    return true;
}

// Do what is necessary if the mesh has moved
void Foam::pointVolInterpolation::topoChange(const polyTopoChangeMap&)
{
    clearAddressing();
    clearGeom();
}

void Foam::pointVolInterpolation::mapMesh(const polyMeshMap&)
{
    clearAddressing();
    clearGeom();
}

void Foam::pointVolInterpolation::distribute(const polyDistributionMap&)
{
    clearAddressing();
    clearGeom();
}


// ************************************************************************* //

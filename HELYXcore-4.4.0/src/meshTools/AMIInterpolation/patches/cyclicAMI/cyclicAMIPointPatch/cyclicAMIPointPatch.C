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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatch/cyclicAMIPointPatch.H"
#include "meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch
    );
    addNamedToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch,
        cyclicPeriodicAMI
    );
    addNamedToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch,
        discreteMixingPlane
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicAMIPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::cyclicAMIPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::cyclicAMIPointPatch::initTopoChange(PstreamBuffers& pBufs)
{
    facePointPatch::initTopoChange(pBufs);
//    cyclicAMIPointPatch::initCalcGeometry(pBufs);
}


void Foam::cyclicAMIPointPatch::topoChange(PstreamBuffers& pBufs)
{
    facePointPatch::topoChange(pBufs);
//    cyclicAMIPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::cyclicAMIPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    cyclicAMIPolyPatch_(refCast<const cyclicAMIPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::~cyclicAMIPointPatch()
{}


// ************************************************************************* //

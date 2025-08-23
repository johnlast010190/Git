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
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunction.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(motionFunction, 0);
    defineRunTimeSelectionTable(motionFunction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::motionFunction::motionFunction
(
    const fvMesh& mesh,
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    mesh_(mesh),
    DBMFCoeffs_(SBMFCoeffs),
    time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::motionFunction::~motionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::motionFunction::boundaryVelocity
(
    const label& pMaster
) const
{
    FatalErrorInFunction
        << "This function should not be called. "
        << "Define them in the derived classes"
        << exit(FatalError);

    return tmp<vectorField>(new vectorField(0, Zero));
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::motionFunction::interfaceVelocity
(
    const label& pMaster
) const
{
    FatalErrorInFunction
        << "This function should not be called "
        << "Define them in the derived classes"
        << exit(FatalError);

    return tmp<vectorField>(new vectorField(0, Zero));
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::motionFunction::interfacePointsVelocity
(
    const label& pMaster
) const
{
    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(mesh_.boundary()[pMaster].patch());
    const labelList& pPoints = gibPolyPatch.meshPoints();

    indirectPatchInterpolation pInterpolation(gibPolyPatch);

    const vectorField U(interfaceVelocity(pMaster));

    // Point interpolated values
    vectorField pointsU(pInterpolation.faceToPointInterpolate(U));
    const scalarField pointsMagU
    (
        pInterpolation.faceToPointInterpolate(mag(U))
    );

    pointsU = pointsMagU*pointsU/(mag(pointsU) + SMALL);

    const vectorField& points = mesh_.points();
    tmp<vectorField> points0(new vectorField(points));

    forAll(pPoints, pI)
    {
        const label pointi = pPoints[pI];
        points0.ref()[pointi] =
            points[pointi] - pointsU[pI]*mesh_.time().deltaTValue();
    }

    return points0;
}


// ************************************************************************* //

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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/morphModesMotion/morphModesMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "meshes/pointMesh/pointMesh.H"
#include "referenceFrames/dynamicMotionCoordinateFrame/dynamicMotionCoordinateFrame.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(morphModesMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        morphModesMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::morphModesMotion::morphModesMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    mesh_(fvSolutionRegistry::getMesh(obr)),
    fileNames_(SBMFCoeffs_.lookup<wordReList>("files")),
    amplifications_(fileNames_.size()),
    triangleModes_(3, vectorField(fileNames_.size(), Zero)),
    CofR_(Zero),
    CofR0_(Zero),
    triPoints_(3, Zero),
    triPoints0_(3, Zero)
{
    morphModesMotion::read(SBMFCoeffs);
    frame().resetDynamic(true);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::morphModesMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    septernion TR = septernion::I;
    if (t2 == 0)
    {
        return TR;
    }

    vectorField deltaMove(3, Zero);
    forAll(triangleModes_, pointi)
    {
        forAll(triangleModes_[pointi], mi)
        {
            deltaMove[pointi] +=
                triangleModes_[pointi][mi]*amplifications_[mi].value(t2);
        }
    }

    // For no motion
    if (max(mag(deltaMove)) < SMALL)
    {
        return TR;
    }

    CofR_ = isIncrementalMotion() ? triPoints_[0] : triPoints0_[0];

    // Calculate two planes normals and angle/axis between them to compute
    // quaternion rotation.
    const vector v1(triPoints_[1] - triPoints_[0]);
    const vector v2(triPoints_[2] - triPoints_[0]);
    const vector n1(normalised(v1^v2));

    // Move here has to be consistent with the move of origin of the frame.
    // That should be ensured by returning correct transformation septernion
    vectorField newPoints(triPoints_ + deltaMove);

    const vector vMoved1(newPoints[1] - newPoints[0]);
    const vector vMoved2(newPoints[2] - newPoints[0]);
    const vector n2(normalised(vMoved1^vMoved2));

    // Rotation around axis (vector line intersection of two planes)
    const vector axis(normalised(n1^n2));

    // No rotation only translating
    if (mag(axis) < SMALL)
    {
        quaternion R(1);

        // When there isn't any rotation the translation for all the points
        // is the same so we can use any vestor from deltaMove.
        TR = septernion(septernion(-deltaMove.first())*R);
    }
    else
    {
        const scalar angle = Foam::acos(n1&n2);
        const quaternion R(axis, angle);
        TR = septernion(-newPoints[0])*R*septernion(triPoints_[0]);
    }

    triPoints_ = isIncrementalMotion() ? newPoints : triPoints0_;

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::labelList
Foam::solidBodyMotionFunctions::morphModesMotion::findTrianglePoints
(
    const vectorField& points
)
{
    const scalarField distanceField(mag(points - frame().coorSys().origin()));
    scalarField threeDistances(3, 1e10);

    // Single precision build will translate (-1) to int64_t buit is label is
    // int32_t, so we need to cast it to label.
    labelList pointIds(3, label(-1));

    List<vectorField> triangleModesValues(3);

    // TODO Looping three times through the field seems like clamsy way to do
    // this quite likely there is a way to do this in one loop.
    // Plus no optimisation is done here, so it can be slow on startup.
    // Of course, there are better methods for finding closest points to other
    // point, like voxalization but due to time constrain it hasn't been
    // considered.

    // Pick up closest point to origin
    forAll(distanceField, i)
    {
        if (distanceField[i] <= threeDistances[0])
        {
            triPoints_[0] = points[i];
            threeDistances[0] = distanceField[i];
            pointIds[0] = i;
        }
    }

    // Pick up second closest point to origin
    forAll(distanceField, i)
    {
        if
        (
            distanceField[i] <= threeDistances[1]
         && distanceField[i] > threeDistances[0]
        )
        {
            triPoints_[1] = points[i];
            threeDistances[1] = distanceField[i];
            pointIds[1] = i;
        }
    }

    // Pick up third closest point to origin
    forAll(distanceField, i)
    {
        if
        (
            distanceField[i] <= threeDistances[2]
         && distanceField[i] > threeDistances[1]
        )
        {
            triPoints_[2] = points[i];
            threeDistances[2] = distanceField[i];
            pointIds[2] = i;
        }
    }
    return pointIds;
}


bool Foam::solidBodyMotionFunctions::morphModesMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    fileNames_ = SBMFCoeffs_.lookup<wordReList>("files");

    //- The modes
    PtrList<pointVectorField> modes(fileNames_.size());
    forAll(fileNames_, fi)
    {
        modes.set
        (
            fi,
            new pointVectorField
            (
                IOobject
                (
                    fileNames_[fi],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                pointMesh::New(mesh_)
            )
        );
    }

    // Find the closest face and three points
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    labelHashSet patchIDs =
        pbm.patchSet(SBMFCoeffs_.lookup<wordReList>("patches"));

    // Collect face centers from all patches and values at points from all modes
    vectorField points;
    List<vectorField> modesBoundaryFields(modes.size());
    if (Pstream::parRun())
    {
        // TODO: We could reduce memory peak by doing distance calculation on
        // each processor separatelly and than collecting just three closest
        // points from each processor choosing which are closest to origin
        // globally and scattering result to all processors.
        // (requires additional work)
        List<vectorField> procPointLoc(Pstream::nProcs());
        List<List<vectorField>> procModesFields(Pstream::nProcs());

        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            const label patchi = *iter;
            procPointLoc[Pstream::myProcNo()] =
                mesh_.boundaryMesh()[patchi].localPoints();
            forAll(modes, mi)
            {
                procModesFields[Pstream::myProcNo()].setSize(modes.size());
                procModesFields[Pstream::myProcNo()][mi] =
                    modes[mi].boundaryField()[patchi].patchInternalField();
            }
            Pstream::allGatherList(procPointLoc);
            Pstream::allGatherList(procModesFields);
            points.append
            (
                ListListOps::combine<vectorField>
                (
                    procPointLoc,
                    accessOp<vectorField>()
                )
            );

            forAll(modes, mi)
            {
                forAll(procModesFields, proci)
                {
                    modesBoundaryFields[mi].append(procModesFields[proci][mi]);
                }
            }
        }
    }
    else
    {
        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            const label patchi = *iter;
            points.append(mesh_.boundaryMesh()[patchi].localPoints());
            forAll(modes, mi)
            {
                modesBoundaryFields[mi].append
                (
                    modes[mi].boundaryField()[patchi].patchInternalField()
                );
            }
        }
    }

    const labelList pointIds = findTrianglePoints(points);
    forAll(triangleModes_, i)
    {
        forAll(modes, mi)
        {
            triangleModes_[i][mi] = modesBoundaryFields[mi][pointIds[i]];
        }
    }
    triPoints0_ = triPoints_;

    CofR_ = triPoints0_[0];
    CofR0_ = triPoints0_[0];

    forAll(fileNames_, fi)
    {
        amplifications_.set
        (
            fi,
            Function1<scalar>::New(fileNames_[fi], SBMFCoeffs_)
        );
    }

    return true;
}


// ************************************************************************* //

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
    (c) 2025 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMotionSolvers/displacement/deformingBodyMotion/deformingBodyFvMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "pointPatchDist/pointPatchDist.H"
#include "interpolation/volPointInterpolation/pointConstraints.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "primitives/functions/Function1/One/OneConstant.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "interpolation/pointVolInterpolation/pointVolInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformingBodyFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        deformingBodyFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        deformingBodyFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::scalar>& Foam::deformingBodyFvMotionSolver::weights
(
    const label pointi,
    List<scalar>& w
) const
{
    // Initialise to 1 for the far-field weight
    scalar sum1mw = 1;

    forAll(bodyMeshes_, bi)
    {
        w[bi] = bodyMeshes_[bi].weight_[pointi];
        sum1mw += w[bi]/(1 + SMALL - w[bi]);
    }

    // Calculate the limiter for wi/(1 - wi) to ensure the sum(wi) = 1
    scalar lambda = 1/sum1mw;

    // Limit wi/(1 - wi) and sum the resulting wi
    scalar sumw = 0;
    forAll(bodyMeshes_, bi)
    {
        w[bi] = lambda*w[bi]/(1 + SMALL - w[bi]);
        sumw += w[bi];
    }

    // Calculate the weight for the stationary far-field
    w[bodyMeshes_.size()] = 1 - sumw;

    return w;
}


void Foam::deformingBodyFvMotionSolver::initBodies
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const dictionary& bodiesDict = dict.subDict("bodies");

    forAllConstIter(IDLList<entry>, bodiesDict, iter)
    {
        const dictionary& bodyDict = iter().dict();

        if (bodyDict.found("patches"))
        {
            bodyMeshes_.append
            (
                new bodyMesh
                (
                    mesh,
                    iter().keyword(),
                    bodyDict
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Body " << iter().keyword()
                << " has no entry patches"
                << exit(FatalError);
        }
    }

    const pointMesh& pMesh = pointMesh::New(mesh);

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        const pointPatchDist pDist(pMesh, bodyMeshes_[bi].patchSet_, points0());

        bodyMeshes_[bi].weight_.primitiveFieldRef() =
            bodyMeshes_[bi].weight(pDist.primitiveField());

        pointConstraints::New(pMesh).constrain(bodyMeshes_[bi].weight_);
    }
}


void Foam::deformingBodyFvMotionSolver::initFields()
{
    if (writeFields_)
    {
        cellDisplacement_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "cellDisplacement",
                    fvMesh_.time().timeName(),
                    fvMesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvMesh_,
                dimensionedVector
                (
                    pointDisplacement_.dimensions(),
                    Zero
                ),
                calculatedFvPatchField<vector>::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deformingBodyFvMotionSolver::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(dict.lookup<scalar>("innerDistance")),
    do_(dict.lookup<scalar>("outerDistance")),
    weight_
    (
        IOobject
        (
            name_ + ".motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, 0)
    )
{}


Foam::deformingBodyFvMotionSolver::deformingBodyFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    writeFields_(coeffDict().lookupOrDefault<Switch>("writeFields", false))
{
    initBodies(mesh,coeffDict());
    initFields();
}


Foam::deformingBodyFvMotionSolver::
deformingBodyFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointVectorField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    writeFields_(coeffDict().lookupOrDefault<Switch>("writeFields", false))
{
    initBodies(mesh,coeffDict());
    initFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::deformingBodyFvMotionSolver::
~deformingBodyFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::deformingBodyFvMotionSolver::bodyMesh::weight
(
    const Type& pDist
) const
{
    // Scaling: 1 up to di then linear down to 0 at do away from patches
    Type weight(min(max((do_ - pDist)/(do_ - di_), scalar(0)), scalar(1)));

    // Convert the weight function to a cosine
    weight =
        min
        (
            max
            (
                0.5 - 0.5*cos(weight*Foam::constant::mathematical::pi),
                scalar(0)
            ),
            scalar(1)
        );

    return weight;
}


Foam::tmp<Foam::pointField>
Foam::deformingBodyFvMotionSolver::curPoints() const
{
    return points0() + pointDisplacement_.primitiveField();
}


Foam::tmp<Foam::pointField>
Foam::deformingBodyFvMotionSolver::analyticalFieldDisplacement
(
    const bodyMesh& body
)
{
    // create pointDisplacement field for body N
    tmp<pointField> tpointDisplacement
    (
        new pointField(pointDisplacement_.primitiveField())
    );
    pointField& pointDisplacement = tpointDisplacement.ref();

    // move all points based on patch motion descriptor
    // (needs bc buffer size to be > domain size)
    forAllConstIter(labelHashSet, body.patchSet_, iter)
    {
        label patchi = iter.key();
        pointDisplacement_.boundaryFieldRef()[patchi].setField
        (
            points0(),
            pointDisplacement,
            true
        );
        // assume the motion for each part of the body is the same
        continue;
    }

    return tpointDisplacement;
}


void Foam::deformingBodyFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    // Update pointDisplacement boundary to new position
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    vectorField& pointDisplacement = pointDisplacement_.primitiveFieldRef();
    const pointField& points0 = this->points0();

    // Update the displacement of internal points
    if (bodyMeshes_.size() == 1)
    {
        const scalarField& weight = bodyMeshes_[0].weight_;
        const pointField transform0
        (
            analyticalFieldDisplacement(bodyMeshes_[0])
        );

        // Weighted interpolation (identity trafo = zero displacement)
        pointDisplacement = transform0 * weight;
    }
    else
    {
        List<scalar> w(bodyMeshes_.size()+1);
        List<pointField> transforms0(bodyMeshes_.size(), pointDisplacement);

        forAll(bodyMeshes_, bi)
        {
            transforms0[bi] = analyticalFieldDisplacement(bodyMeshes_[bi]);
        }

        forAll(points0, pointi)
        {
            // update weights for point i
            weights(pointi, w);

            // get displacements for all bodies (N+1 displacement = zero)
            List<point> transform0(bodyMeshes_.size()+1, vector::zero);
            forAll(bodyMeshes_, bi)
            {
                transform0[bi] = transforms0[bi][pointi];
            }

            // average operator does not exist for Lists
            // add displacement contributions of each body manually
            pointDisplacement[pointi] = transform0[0]*w[0];
            for (label i=1; i<transform0.size(); i++)
            {
                pointDisplacement[pointi] += transform0[i]*w[i];
            }
        }
    }

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);

    // Optional compute cellDisplacement
    if (writeFields_)
    {
        pointVolInterpolation::New(fvMesh_).interpolate
        (
            pointDisplacement_,
            cellDisplacement_()
        );
    }
}


void Foam::deformingBodyFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    displacementMotionSolver::topoChange(map);
}


void Foam::deformingBodyFvMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    displacementMotionSolver::mapMesh(map);
}


// ************************************************************************* //

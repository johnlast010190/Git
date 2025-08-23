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
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/leastSquaresGrad/leastSquaresVectors.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcCellReduce.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresVectors>(mesh),
    pVectorsPtr_(nullptr),
    nVectorsPtr_(nullptr),
    wLayerCellsPtr_(nullptr)
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
    deleteDemandDrivenData(wLayerCellsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectors::calcLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresVectors::calcLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, Zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, Zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const surfaceScalarField& w = mesh_.weights();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), Zero);

    tmp<surfaceVectorField> td = mesh_.delta();
    const surfaceVectorField& d = td();
    tmp<surfaceVectorField> tmd = fvc::applyFaceMask(d);
    const surfaceVectorField& md = tmd();
    surfaceSymmTensorField mdSqr(sqr(d));
    fvc::applyFaceMaskTo(mdSqr);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        symmTensor wdd = (1.0/magSqr(d[facei]))*mdSqr[facei];

        dd[own] += wdd;
        dd[nei] += wdd;
    }


    surfaceVectorField::Boundary& pVectorsBf =
        lsP.boundaryFieldRef();

    const volVectorField& C = mesh_.C();

    forAll(pVectorsBf, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        if (!pw.coupled())
        {
            const vectorField& pcf = pw.patch().Cf();
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                vector dist = pcf[patchFacei]-C[own];
                dd[own] += (1.0/magSqr(dist))*sqr(dist);
            }
        }
        else
        {
            const vectorField& pd = d.boundaryField()[patchi];
            const symmTensorField& pmdSqr = mdSqr.boundaryField()[patchi];

            forAll(pd, patchFacei)
            {
                dd[faceCells[patchFacei]] +=
                    (1.0/magSqr(pd[patchFacei]))*pmdSqr[patchFacei];
            }
        }
    }


    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));

    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar magSfByMagSqrd = 1.0/magSqr(d[facei]);

        lsP[facei] = magSfByMagSqrd*(invDd[own] & md[facei]);
        lsN[facei] = -magSfByMagSqrd*(invDd[nei] & md[facei]);
    }

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();
        const vectorField& pcf = pw.patch().Cf();

        if (!pw.coupled())
        {
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                vector dist = pcf[patchFacei]-C[own];
                patchLsP[patchFacei] = (1.0/magSqr(dist))
                    *(invDd[own] & dist);
            }
        }
        else
        {
            const vectorField& pd = d.boundaryField()[patchi];
            const vectorField& pmd = md.boundaryField()[patchi];

            forAll(pd, patchFacei)
            {
                patchLsP[patchFacei] =
                    (1.0/magSqr(pd[patchFacei]))
                    *(invDd[faceCells[patchFacei]] & pmd[patchFacei]);
            }
        }
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing least square gradient vectors" << endl;
    }
}


void Foam::leastSquaresVectors::calcWallLayerCells(const label nLayers) const
{
    if (debug)
    {
        Info<< "leastSquaresVectors::calcWallLayerCells() :"
            << "Calculating layer cell field indentifier"
            << endl;
    }

    wLayerCellsPtr_ = new volScalarField
    (
        IOobject
        (
            "LeastSquaresLayerCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero),
        "zeroGradient"
    );
    volScalarField& wC = *wLayerCellsPtr_;
    forAll(this->mesh().boundary(), pI)
    {
        if (isA<wallFvPatch>(this->mesh().boundary()[pI]))
        {
            const labelUList& faceCells
            (
                this->mesh().boundary()[pI].faceCells()
            );

            forAll(faceCells, fcI)
            {
                wC.operator[](faceCells[fcI]) = 1;
            }
        }
    }
    wC.correctBoundaryConditions();

    for (label ai = 2; ai <= nLayers; ai++)
    {
        string localMax("localMax");
        surfaceScalarField layerf
        (
            fvc::interpolate(wC, IStringStream(localMax)())
        );
        wC = fvc::cellReduce(layerf, maxEqOp<scalar>(), scalar(0.0));
        wC.correctBoundaryConditions();
    }
    if (false)
    {
        wC.write();
    }

    if (debug)
    {
        InfoInFunction
            << "Finished calculating layer cells" << endl;
    }
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        calcLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        calcLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}

const Foam::volScalarField& Foam::leastSquaresVectors::wallLayerCells
(
    const label layers
)
const
{
    if (!wLayerCellsPtr_)
    {
        calcWallLayerCells(layers);
    }

    return *wLayerCellsPtr_;
}


bool Foam::leastSquaresVectors::movePoints()
{
    if (debug)
    {
        InfoIn("bool leastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
    deleteDemandDrivenData(wLayerCellsPtr_);

    return true;
}



// ************************************************************************* //

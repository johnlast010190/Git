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
    (c) 2011 OpenFOAM Foundation
    (c) 2016-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/taylorGauss/taylorGaussData.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(taylorGaussData, 0);

float taylorGaussData::condNuThreshold_
(
    Foam::debug::floatOptimisationSwitch("taylorGaussThreshold", 10)
);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::tensorField> Foam::taylorGaussData::pseudoVolumes() const
{
    if (debug)
    {
        Info<< "taylorGaussNormals::makeInvPseudoVolumes() :"
            << "Constructing inverse pseudo-volume tensors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();

    // Set up temporary storage for edge centres and calculate
    surfaceVectorField Ce
    (
        IOobject
        (
            "Ce",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimLength
    );

    tmp<tensorField> pseudoVolPtr
    (
        new tensorField(mesh_.nCells(), tensor::zero)
    );

    tensorField& pseudoVol(pseudoVolPtr.ref());

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
    const surfaceVectorField& mSf = tmSf();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        Ce[facei] = w[facei]*C[own] + (1-w[facei])*C[nei];
        pseudoVol[own] += mSf[facei]*(Ce[facei]-C[own]);
        pseudoVol[nei] -= mSf[facei]*(Ce[facei]-C[nei]);
    }

    forAll(w.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        fvsPatchVectorField& pCe = Ce.boundaryFieldRef()[patchi];
        const fvsPatchVectorField& pmSf = mSf.boundaryField()[patchi];
        const labelUList& faceCells = pw.patch().faceCells();

        if (pw.coupled())
        {
            tmp<vectorField> tpd = pw.patch().delta();
            const vectorField& pd = tpd();

            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                pCe[patchFacei] = C[own]+(1-pw[patchFacei])*pd[patchFacei];
                pseudoVol[own] += pmSf[patchFacei]*(pCe[patchFacei]-C[own]);
            }
        }
        else
        {
            const vectorField& pcf = pw.patch().Cf();
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                pseudoVol[own] += pmSf[patchFacei]*(pcf[patchFacei]-C[own]);
            }
        }
    }

    mesh.stabiliseEmptyDirections(pseudoVol);

    if (debug)
    {
        Info<< "taylorGaussNormals::makeCoeffs() :"
            << "Finished constructing skew Gauss data"
            << endl;
    }

    return pseudoVolPtr;
}


template<class Type>
Foam::tmp<Foam::tensorField> Foam::taylorGaussData::pseudoVolumes
(
    VolField<Type>  vf
) const
{
    const fvMesh& mesh(vf.mesh());

    // treatment for extrapolated boundaries
    tmp<tensorField> pseudoVolPtr(pseudoVolumes());
    tensorField& pseudoVol(pseudoVolPtr.ref());

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh_.Sf());
    const surfaceVectorField& mSf = tmSf();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);

        const labelUList& pFaceCells = patch.faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        if (vf.boundaryField()[patchi].extrapolated())
        {
            tmp<vectorField> delta(patch.delta());

            forAll(patch, facei)
            {
                pseudoVol[pFaceCells[facei]]
                    -= (pmSf[facei]*delta->operator[](facei));
            }
        }
    }

    mesh.stabiliseEmptyDirections(pseudoVol);

    return pseudoVolPtr;
}


Foam::tmp<Foam::scalarField> Foam::taylorGaussData::calcConditionNumber() const
{
    if (debug)
    {
        Info<< "taylorGaussNormals::calcConditionNumber() :"
            << "Constructing condition number field of pseudoVolume tensors"
            << endl;
    }

    tmp<scalarField> condNumPtr
    (
        new scalarField(mesh_.nCells(), 1.0)
    );

    scalarField& condNum(condNumPtr.ref());

    //- pseudoVol call. Called two times/moveMesh. One for inv one for condNu
    //  Need refactoring in the future
    tmp<tensorField> tpseudoVol = pseudoVolumes()/mesh().V();
    tmp<symmTensorField> tpseudoVol2
        (
            symm
            (
                 (tpseudoVol() & tpseudoVol().T())
            )
        );
    symmTensorField& pseudoVol2 = tpseudoVol2.ref();
    //- crop small values to make eigenvalue algorithm better behaved
    forAll(pseudoVol2, cI)
    {
        if (mag(pseudoVol2[cI].xy()) < 1e-6) pseudoVol2[cI].xy() = 0;
        if (mag(pseudoVol2[cI].xz()) < 1e-6) pseudoVol2[cI].xz() = 0;
        if (mag(pseudoVol2[cI].yz()) < 1e-6) pseudoVol2[cI].yz() = 0;
        if (mag(pseudoVol2[cI].xx()) < 0.0)  pseudoVol2[cI].xx() = 0;
        if (mag(pseudoVol2[cI].yy()) < 0.0)  pseudoVol2[cI].yy() = 0;
        if (mag(pseudoVol2[cI].zz()) < 0.0)  pseudoVol2[cI].zz() = 0;
    }

    label nSol = mesh_.nSolutionD();
    label cpntMax = 2;
    if (nSol==2) cpntMax = 1;


    //- calculate singular values (sigma) of pseudoVol and then
    //  sigmaMax/sigmaMin
    if (nSol>1)
    {
        forAll(condNum, cI)
        {
            vector eigens = eigenValues(pseudoVol2[cI]);
            //- Usually they return sorted - check tensor.H
            scalar meigenMin = eigens.component(0);
            scalar meigenMax = eigens.component(cpntMax);
            if ((meigenMax-meigenMin)>SMALL)
            {
                if (meigenMin > ROOTVSMALL)
                {
                    condNum[cI] = meigenMax/meigenMin;
                }
            }
        }
    }

    return condNumPtr;
}


void Foam::taylorGaussData::clearData()
{
    deleteDemandDrivenData(invPseudoVolPtr_);
    deleteDemandDrivenData(fieldInvPseudoVolPtrs_);
    deleteDemandDrivenData(cAPtr_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::taylorGaussData::taylorGaussData(const fvMesh& mesh)
:
    MeshObject<fvMesh, MoveableMeshObject, taylorGaussData>(mesh),
    invPseudoVolPtr_(nullptr),
    fieldInvPseudoVolPtrs_
    (
        new HashPtrTable<Field<tensor>, word, string::hash>()
    ),
    cAPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::taylorGaussData::~taylorGaussData()
{
    clearData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes() const
{
    if (!invPseudoVolPtr_)
    {
        // Store inverse of tensor
        invPseudoVolPtr_ = inv(pseudoVolumes()).ptr();
    }

    return *invPseudoVolPtr_;
}


const Foam::scalarField& Foam::taylorGaussData::conditionNumber() const
{
    if (!cAPtr_)
    {
        // Calculate and store condition number
        cAPtr_ = calcConditionNumber().ptr();
    }

    return *cAPtr_;
}


const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes
(
    const volScalarField& vsf
) const
{
    // check if field specific modification as required
    bool extrapolated = false;
    forAll(vsf.boundaryField(), patchi)
    {
        if (vsf.boundaryField()[patchi].extrapolated())
        {
            extrapolated = true;
            break;
        }
    }

    if (extrapolated)
    {
        if (!fieldInvPseudoVolPtrs_->found(vsf.name()))
        {
            fieldInvPseudoVolPtrs_->insert
            (
                vsf.name(),
                (inv(pseudoVolumes(vsf))).ptr()
            );
        }

        return *(fieldInvPseudoVolPtrs_->operator[](vsf.name()));
    }

    return invPseudoVolumes();
}


const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes
(
    const volVectorField& vvf
) const
{
    // check if field specific modification as required
    bool extrapolated = false;
    forAll(vvf.boundaryField(), patchi)
    {
        if (vvf.boundaryField()[patchi].extrapolated())
        {
            extrapolated = true;
            break;
        }
    }

    if (extrapolated)
    {
        if (!fieldInvPseudoVolPtrs_->found(vvf.name()))
        {
            fieldInvPseudoVolPtrs_->insert
            (
                vvf.name(),
                (inv(pseudoVolumes(vvf))).ptr()
            );
        }

        return *(fieldInvPseudoVolPtrs_->operator[](vvf.name()));
    }

    return invPseudoVolumes();
}


bool Foam::taylorGaussData::movePoints()
{
    clearData();

    return true;
}


// ************************************************************************* //

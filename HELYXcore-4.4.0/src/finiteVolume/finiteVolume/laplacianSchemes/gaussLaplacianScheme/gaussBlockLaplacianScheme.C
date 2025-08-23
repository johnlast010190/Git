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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline tmp<fvBlockMatrix<vector>>
gaussLaplacianScheme<vector, tensor>::fvmBLaplacian
(
    const SurfaceField<tensor>& gamma,
    const VolField<vector>& vf
)
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());

    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const SurfaceField<scalar> gammaMagSf
    (
        SfGamma & Sn
    );

    tmp<surfaceScalarField> tdeltaCoeffs
    (
        this->tsnGradScheme_().deltaCoeffs(vf)
    );
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    tmp<fvBlockMatrix<vector>> tfvbm
    (
        new fvBlockMatrix<vector>
        (
            const_cast<VolField<vector>&>(vf),
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );

    fvBlockMatrix<vector>& fvbm = tfvbm.ref();
    vectorField& source = fvbm.source();

    CoeffField<vector>::squareTypeField& d = fvbm.diag().asSquare();
    CoeffField<vector>::squareTypeField& u = fvbm.upper().asSquare();
    CoeffField<vector>::squareTypeField& l = fvbm.lower().asSquare();

    const direction nCmpts = pTraits<vector>::nComponents;
    const scalarField dgMsf
    (
        deltaCoeffs.primitiveField() * gammaMagSf.internalField()
    );

    forAll(u, fI)
    {
        for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            //- lower not needed. Need to make it asymetrix to be able to
            //  run AMG for testing purposes
            u[fI](cmptI, cmptI) = dgMsf[fI];
            l[fI](cmptI, cmptI) = dgMsf[fI];
        }
    }
    fvbm.negSumDiag();

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchField<vector>& pvf = vf.boundaryField()[pI];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[pI];

        const fvPatchVectorField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const labelUList& fc = patch.faceCells();

        if (pvf.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                fvbm.coupleUpper()[pI].asLinear();

            CoeffField<vector>::linearTypeField& pcLower =
                fvbm.coupleLower()[pI].asLinear();

            const fvsPatchScalarField& pDeltaCoeffs =
                deltaCoeffs.boundaryField()[pI];

            //- For now we can use only the segregated coeffs
            //  In the future for AMI etc this has to be changed
            const vectorField iCoeffs(pvf.gradientInternalCoeffs(pDeltaCoeffs));
            const vectorField bCoeffs(pvf.gradientBoundaryCoeffs(pDeltaCoeffs));

            forAll(pf, fI)
            {
                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    d[fc[fI]](cmptI, cmptI) += pGamma[fI]*iCoeffs[fI](cmptI);
                }
            }

            pcUpper -= pGamma*bCoeffs;
            pcLower += pGamma*iCoeffs;
        }
        else
        {
            CoeffField<vector> biCoeffs(pvf.gradientInternalBCoeffs());
            const CoeffField<vector>::squareTypeField& sbiCoeffs =
                biCoeffs.asSquare();

            const vectorField bbCoeffs(pvf.gradientBoundaryBCoeffs());

            forAll(pf, fI)
            {
                const label fcI = fc[fI];

                d[fcI] += pGamma[fI]*sbiCoeffs[fI];
                source[fcI] -= pGamma[fI]*bbCoeffs[fI];
            }
        }
    }

    if (this->tsnGradScheme_().corrected())
    {
        if (mesh.schemes().fluxRequired(vf.name()))
        {
            surfaceVectorField faceFluxCorr ( gammaMagSf*this->tsnGradScheme_().correction(vf) );

            source -= mesh.V()*
                fvc::div(faceFluxCorr)().primitiveField();
        }
        else
        {
            source -= mesh.V()*
                fvc::div
                (
                    gammaMagSf*this->tsnGradScheme_().correction(vf)
                )().primitiveField();
        }
    }



    return tfvbm;
}


template<class Type, class GType>
inline tmp<fvBlockMatrix<Type>>
gaussLaplacianScheme<Type, GType>::fvmBLaplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf
)
{
    FatalErrorInFunction
        << "BLaplacianscalar(scalar, scalar-vector) suppoted"
        << exit(FatalError);

    tmp<fvBlockMatrix<Type>> tfvbm
    (
        new fvBlockMatrix<Type>
        (
            const_cast<VolField<Type>&>(vf)
        )
    );

    return tfvbm;
}


template<>
inline tmp<fvBlockMatrix<scalar>>
gaussLaplacianScheme<scalar, scalar>::fvmBLaplacian
(
    const SurfaceField<scalar>& gamma,
    const VolField<scalar>& vf
)
{
    const fvMesh& mesh = this->mesh();

    SurfaceField<scalar> gammaMagSf
    (
        gamma*mesh.magSf()
    );

    tmp<surfaceScalarField> tdeltaCoeffs
    (
        this->tsnGradScheme_().deltaCoeffs(vf)
    );
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    tmp<fvBlockMatrix<scalar>> tfvbm
    (
        new fvBlockMatrix<scalar>
        (
            const_cast<VolField<scalar>&>(vf)
        )
    );
    fvBlockMatrix<scalar>& fvbm = tfvbm.ref();
    scalarField& source = fvbm.source();

    CoeffField<scalar>::linearTypeField& d = fvbm.diag().asScalar();
    CoeffField<scalar>::linearTypeField& u = fvbm.upper().asScalar();

    u = deltaCoeffs.primitiveField()*gammaMagSf.internalField();
    fvbm.negSumDiag();

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchField<scalar>& pvf = vf.boundaryField()[pI];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[pI];

        const fvPatchScalarField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const labelUList& fc = patch.faceCells();

        if (pvf.coupled())
        {
            CoeffField<scalar>::linearTypeField& pcUpper =
                fvbm.coupleUpper()[pI].asScalar();

            CoeffField<scalar>::linearTypeField& pcLower =
                fvbm.coupleLower()[pI].asLinear();

            const fvsPatchScalarField& pDeltaCoeffs =
                deltaCoeffs.boundaryField()[pI];

            const scalarField iCoeffs(pvf.gradientInternalCoeffs(pDeltaCoeffs));
            const scalarField bCoeffs(pvf.gradientBoundaryCoeffs(pDeltaCoeffs));

            forAll(pf, fI)
            {
                d[fc[fI]] += pGamma[fI]*iCoeffs[fI];
            }

            pcUpper -= pGamma*bCoeffs;
            pcLower += pGamma*iCoeffs;
        }
        else
        {
            const scalarField iCoeffs(pvf.gradientInternalCoeffs());
            const scalarField bCoeffs(pvf.gradientBoundaryCoeffs());

            forAll(pf, fI)
            {
                const label fcI = fc[fI];

                d[fcI] += pGamma[fI]*iCoeffs[fI];
                source[fcI] -= pGamma[fI]*bCoeffs[fI];
            }
        }
    }

    if (this->tsnGradScheme_().corrected())
    {
        if (mesh.schemes().fluxRequired(vf.name()))
        {
            surfaceScalarField faceFluxCorr
            (
                gammaMagSf*this->tsnGradScheme_().correction(vf)
            );

            source -= mesh.V()*
                fvc::div(faceFluxCorr)().primitiveField();
        }
        else
        {
            source -= mesh.V()*
                fvc::div
                (
                    gammaMagSf*this->tsnGradScheme_().correction(vf)
                )().primitiveField();
        }
    }
    return tfvbm;
}



template<>
inline tmp<fvBlockMatrix<vector>>
gaussLaplacianScheme<vector, scalar>::fvmBLaplacian
(
    const SurfaceField<scalar>& gamma,
    const VolField<vector>& vf
)
{
    const fvMesh& mesh = this->mesh();

    SurfaceField<scalar> gammaMagSf
    (
        gamma*mesh.magSf()
    );

    tmp<surfaceScalarField> tdeltaCoeffs
    (
        this->tsnGradScheme_().deltaCoeffs(vf)
    );
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    tmp<fvBlockMatrix<vector>> tfvbm
    (
        new fvBlockMatrix<vector>
        (
            const_cast<VolField<vector>&>(vf),
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );

    tmp<surfaceScalarField> tmGammaMagSf = fvc::applyFaceMask(gammaMagSf);
    const surfaceScalarField& mGammaMagSf = tmGammaMagSf();

    fvBlockMatrix<vector>& fvbm = tfvbm.ref();
    vectorField& source = fvbm.source();

    CoeffField<vector>::squareTypeField& d = fvbm.diag().asSquare();
    CoeffField<vector>::squareTypeField& u = fvbm.upper().asSquare();
    CoeffField<vector>::squareTypeField& l = fvbm.lower().asSquare();

    const direction nCmpts = pTraits<vector>::nComponents;
    const scalarField dgMsf
    (
        deltaCoeffs.primitiveField() * mGammaMagSf.internalField()
    );

    forAll(u, fI)
    {
        for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            //- lower not needed. Need to make it asymetrix to be able to
            //  run AMG for testing purposes
            u[fI](cmptI, cmptI) = dgMsf[fI];
            l[fI](cmptI, cmptI) = dgMsf[fI];
        }
    }
    fvbm.negSumDiag();

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchField<vector>& pvf = vf.boundaryField()[pI];
        const fvsPatchScalarField& pGamma = mGammaMagSf.boundaryField()[pI];

        const fvPatchVectorField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const labelUList& fc = patch.faceCells();

        if (pvf.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                fvbm.coupleUpper()[pI].asLinear();

            CoeffField<vector>::linearTypeField& pcLower =
                fvbm.coupleLower()[pI].asLinear();

            const fvsPatchScalarField& pDeltaCoeffs =
                deltaCoeffs.boundaryField()[pI];

            //- For now we can use only the segregated coeffs
            //  In the future for AMI etc this has to be changed
            const vectorField iCoeffs(pvf.gradientInternalCoeffs(pDeltaCoeffs));
            const vectorField bCoeffs(pvf.gradientBoundaryCoeffs(pDeltaCoeffs));

            forAll(pf, fI)
            {
                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    d[fc[fI]](cmptI, cmptI) += pGamma[fI]*iCoeffs[fI](cmptI);
                }
            }

            pcUpper -= pGamma*bCoeffs;
            pcLower += pGamma*iCoeffs;
        }
        else
        {
            CoeffField<vector> biCoeffs(pvf.gradientInternalBCoeffs());
            const CoeffField<vector>::squareTypeField& sbiCoeffs =
                biCoeffs.asSquare();

            const vectorField bbCoeffs(pvf.gradientBoundaryBCoeffs());

            forAll(pf, fI)
            {
                const label fcI = fc[fI];

                d[fcI] += pGamma[fI]*sbiCoeffs[fI];
                source[fcI] -= pGamma[fI]*bbCoeffs[fI];
            }
        }
    }

    if (this->tsnGradScheme_().corrected())
    {
        if (mesh.schemes().fluxRequired(vf.name()))
        {
            surfaceVectorField faceFluxCorr
            (
                mGammaMagSf*this->tsnGradScheme_().correction(vf)
            );

            source -= mesh.V()*
                fvc::div(faceFluxCorr)().primitiveField();
        }
        else
        {
            source -= mesh.V()*
                fvc::div
                (
                    mGammaMagSf*this->tsnGradScheme_().correction(vf)
                )().primitiveField();
        }
    }


    return tfvbm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(gaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvmLaplacian         \
(                                                                              \
    const SurfaceField<scalar>& gamma,                                         \
    const VolField<Type>& vf                                                   \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    SurfaceField<scalar> gammaMagSf                                            \
    (                                                                          \
        gamma*mesh.magSf()                                                     \
    );                                                                         \
                                                                               \
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected                         \
    (                                                                          \
        gammaMagSf,                                                            \
        this->tsnGradScheme_().deltaCoeffs(vf),                                \
        vf                                                                     \
    );                                                                         \
    fvMatrix<Type>& fvm = tfvm.ref();                                          \
                                                                               \
    if (this->tsnGradScheme_().corrected())                                    \
    {                                                                          \
        if (mesh.schemes().fluxRequired(vf.name()))                            \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = new                                  \
            SurfaceField<Type>                                                 \
            (                                                                  \
                gammaMagSf*this->tsnGradScheme_().correction(vf)               \
            );                                                                 \
                                                                               \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    *fvm.faceFluxCorrectionPtr()                               \
                )().primitiveField();                                          \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)           \
                )().primitiveField();                                          \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::VolField<Foam::Type>> \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian         \
(                                                                              \
    const SurfaceField<scalar>& gamma,                                         \
    const VolField<Type>& vf                                                   \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    tmp<VolField<Type>> tLaplacian                \
    (                                                                          \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())         \
    );                                                                         \
                                                                               \
    tLaplacian.ref().rename                                                    \
    (                                                                          \
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'                    \
    );                                                                         \
                                                                               \
    return tLaplacian;                                                         \
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


#define declareFvmScalarBLaplacianGamma(Type)                                  \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvBlockMatrix<Foam::scalar>>                                   \
Foam::fv::gaussLaplacianScheme<Foam::scalar, Foam::Type>::fvmBLaplacian        \
(                                                                              \
    const SurfaceField<Type>& gamma,             \
    const VolField<scalar>& vf                    \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());                       \
                                                                               \
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);                       \
    const SurfaceField<scalar> gammaMagSf        \
    (                                                                          \
        SfGamma & Sn                                                           \
    );                                                                         \
                                                                               \
    tmp<surfaceScalarField> tdeltaCoeffs                                       \
    (                                                                          \
        this->tsnGradScheme_().deltaCoeffs(vf)                                 \
    );                                                                         \
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();                    \
                                                                               \
    tmp<fvBlockMatrix<scalar>> tfvbm                                           \
    (                                                                          \
        new fvBlockMatrix<scalar>                                              \
        (                                                                      \
            const_cast<VolField<scalar>&>(vf)     \
        )                                                                      \
    );                                                                         \
    fvBlockMatrix<scalar>& fvbm = tfvbm.ref();                                 \
    scalarField& source = fvbm.source();                                       \
                                                                               \
    CoeffField<scalar>::linearTypeField& d = fvbm.diag().asScalar();           \
    CoeffField<scalar>::linearTypeField& u = fvbm.upper().asScalar();          \
                                                                               \
    u = deltaCoeffs.primitiveField()*gammaMagSf.internalField();               \
    fvbm.negSumDiag();                                                         \
                                                                               \
    forAll(vf.boundaryField(), pI)                                             \
    {                                                                          \
        const fvPatchField<scalar>& pvf = vf.boundaryField()[pI];              \
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[pI];    \
        const fvsPatchScalarField& pDeltaCoeffs =                              \
            deltaCoeffs.boundaryField()[pI];                                   \
                                                                               \
        const fvPatchScalarField& pf = vf.boundaryField()[pI];                 \
        const fvPatch& patch = pf.patch();                                     \
        const labelUList& fc = patch.faceCells();                              \
                                                                               \
        const scalarField iCoeffs(pvf.gradientInternalCoeffs(pDeltaCoeffs));   \
                                                                               \
        forAll(pf, fI)                                                         \
        {                                                                      \
            d[fc[fI]] += pGamma[fI]*iCoeffs[fI];                               \
        }                                                                      \
                                                                               \
        if (pvf.coupled())                                                     \
        {                                                                      \
            CoeffField<scalar>::linearTypeField& pcUpper =                     \
                fvbm.coupleUpper()[pI].asScalar();                             \
                                                                               \
            pcUpper -= pGamma*iCoeffs;                                         \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            const scalarField bCoeffs(pvf.gradientBoundaryCoeffs(pDeltaCoeffs));\
                                                                               \
            forAll(pf, fI)                                                     \
            {                                                                  \
                source[fc[fI]] -= pGamma[fI]*bCoeffs[fI];                      \
            }                                                                  \
                                                                               \
        }                                                                      \
    }                                                                          \
    FatalErrorInFunction                                                       \
        << "scalar not suppoted"                                               \
        << abort(FatalError);                                                  \
                                                                               \
    return tfvbm;                                                              \
}                                                                              \


declareFvmScalarBLaplacianGamma(sphericalTensor);
declareFvmScalarBLaplacianGamma(symmTensor);
declareFvmScalarBLaplacianGamma(tensor);

// ************************************************************************* //

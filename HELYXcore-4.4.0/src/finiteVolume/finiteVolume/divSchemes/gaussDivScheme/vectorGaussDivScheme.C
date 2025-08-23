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

Description
    Specialisation of gaussDivScheme for vectors

Copyright
    (c) 2019-2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/divSchemes/gaussDivScheme/gaussDivScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, scalar>> gaussDivScheme<vector>::fvmDiv
(
    const VolField<vector>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, scalar>> tbs
    (
        new BlockLduSystem<vector, scalar>(mesh)
    );
    BlockLduSystem<vector, scalar>& bs = tbs.ref();
    scalarField& source = bs.source();

    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    surfaceScalarField weights (  this->tinterpScheme_().weights(vf) );

    const surfaceVectorField mSf(fvc::applyFaceMask(mesh.Sf()));

    const vectorField& sf = mSf.internalField();
    const scalarField& w = weights.internalField();

    l = -w*sf;
    u = l + sf;
    bs.negSumDiag();

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchVectorField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = mSf.boundaryField()[pI];
        const fvsPatchScalarField& pw = weights.boundaryField()[pI];
        const labelUList& fc = patch.faceCells();

        const vectorField iCoeffs(pf.valueDivInternalCoeffs(pw));
        const vectorField bCoeffs(pf.valueDivBoundaryCoeffs(pw));

        bs.internalCoeffs().set(pI, new CoeffField<vector>(pf.size()));
        CoeffField<vector>::linearTypeField& iC =
            bs.internalCoeffs()[pI].asLinear();

        iC = iCoeffs;

        bs.boundaryCoeffs().set(pI, new Field<scalar>(pf.size(), 0));
        scalarField& bC = bs.boundaryCoeffs()[pI];

        iC = cmptMultiply(iCoeffs, pSf);

        forAll(pf, fI)
        {
            d[fc[fI]] += iC[fI];
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                bs.coupleUpper()[pI].asLinear();
            CoeffField<vector>::linearTypeField& pcLower =
                bs.coupleLower()[pI].asLinear();

            const vectorField pcl (iC);
            const vectorField pcu (cmptMultiply(bCoeffs, pSf));
            pcLower += pcl;
            pcUpper -= pcu;
        }
        else
        {
            bC = -bCoeffs&pSf;
            forAll(pf, fI)
            {
                source[fc[fI]] += bC[fI];
            }
        }

    }

    if (tinterpScheme_().corrected())
    {
        tmp<surfaceScalarField> tmCorr =
            mSf & tinterpScheme_().correction(vf);
        source -= fvc::surfaceIntegrate(tmCorr())*mesh.V();
        //if (vf.mesh().schemes().fluxRequired(vf.name()))
        // Currently div(U) is always at pressure eq
        {
            bs.faceFluxCorrectionPtr() = tmCorr.ptr();
        }
    }

    return tbs;
}


template<>
tmp<BlockLduSystem<vector, scalar>> gaussDivScheme<vector>::fvmBCDiv
(
    const SurfaceField<scalar>& rhof,
    const VolField<vector>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, scalar>> tbs
    (
        new BlockLduSystem<vector, scalar>(mesh)
    );
    BlockLduSystem<vector, scalar>& bs = tbs.ref();
    scalarField& source = bs.source();

    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    surfaceScalarField weights (  this->tinterpScheme_().weights(vf) );
    surfaceVectorField mrhofSf(fvc::applyFaceMask(rhof*mesh.Sf()));

    const scalarField& w = weights.internalField();
    const vectorField& mrhofSfInt = mrhofSf.internalField();

    l = -w*mrhofSfInt;
    u = l + mrhofSfInt;
    bs.negSumDiag();

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchVectorField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const fvsPatchScalarField& pw = weights.boundaryField()[pI];
        const fvsPatchVectorField& pmRhofSf = mrhofSf.boundaryField()[pI];
        const labelUList& fc = patch.faceCells();

        const vectorField iCoeffs(pf.valueDivInternalCoeffs(pw));
        const vectorField bCoeffs(pf.valueDivBoundaryCoeffs(pw));

        bs.internalCoeffs().set(pI, new CoeffField<vector>(pf.size()));
        CoeffField<vector>::linearTypeField& iC =
            bs.internalCoeffs()[pI].asLinear();

        bs.boundaryCoeffs().set(pI, new Field<scalar>(pf.size(), 0));
        scalarField& bC = bs.boundaryCoeffs()[pI];

        iC = cmptMultiply(iCoeffs, pmRhofSf);

        forAll(pf, fI)
        {
            d[fc[fI]] += iC[fI];
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                bs.coupleUpper()[pI].asLinear();
            CoeffField<vector>::linearTypeField& pcLower =
                bs.coupleLower()[pI].asLinear();

            const vectorField pcl (iC);
            const vectorField pcu (cmptMultiply(bCoeffs, pmRhofSf));
            pcLower += pcl;
            pcUpper -= pcu;
        }
        else
        {
            bC = -bCoeffs&(pmRhofSf);
            forAll(pf, fI)
            {
                source[fc[fI]] += bC[fI];
            }
        }
    }

    if (tinterpScheme_().corrected())
    {
        tmp<SurfaceField<scalar>> tmCorr =
            mrhofSf & tinterpScheme_().correction(vf);
        source -= fvc::surfaceIntegrate(tmCorr())*mesh.V();
        //if (vf.mesh().schemes().fluxRequired(vf.name()))
        // Currently div(U) is always at pressure eq
        {
            bs.faceFluxCorrectionPtr() = tmCorr.ptr();
        }
    }

    return tbs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

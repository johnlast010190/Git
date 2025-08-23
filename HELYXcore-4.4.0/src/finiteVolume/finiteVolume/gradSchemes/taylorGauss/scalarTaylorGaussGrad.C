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

Description
    Specialisation of taylorGaussGrad for scalars

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/taylorGauss/scalarTaylorGaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, vector>> taylorGaussGrad<scalar>::fvmGrad
(
    const VolField<scalar>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, vector>> tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs.ref();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    const vectorField& sf = mesh.Sf().internalField();

    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(vf);
    const scalarField& w = tweights().internalField();

    // Get reference to skew face vectors
    const taylorGaussData& data = taylorGaussData::New(mesh);
    const tensorField& invPseudoVol = data.invPseudoVolumes(vf);
    const scalarField& vol = mesh.V();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    tmp<scalarField> tcondNu = data.conditionNumber();
    const scalarField& condNu = tcondNu();

    boolList markCell(mesh.nCells(), false);
    scalar TGcells(0);

    const scalar condNuThreashold = taylorGaussData::condNuThreshold_;
    forAll(markCell, cI)
    {
        if (condNu[cI] < condNuThreashold)
        {
            markCell[cI] = true;
            TGcells++;
        }
    }
    if (debug && mesh.time().outputTime())
    {
        printConditionNumber(vf.mesh(), condNu);
    }

    if (debug)
    {
        scalar allcells(mesh.nCells());
        reduce(TGcells, sumOp<scalar>());
        reduce(allcells, sumOp<scalar>());
        Info<< "ConditionNumberThreshold: " << condNuThreashold << endl;
        Info<< "TGcells: " << TGcells << endl;
        Info<< "GGcells: " << (allcells - TGcells) << endl;
    }

    //- if defered add TG-GG difference to the RHS
    //  Currently can create instabilities.
    //  Normally the problematic cells are very few, therefore better to be
    //  safe there and avoid the use of inverted quantity.
    bool deferred = false;

    forAll(own, fI)
    {
        const label ownFI = own[fI];
        const label neiFI = nei[fI];
        const vector wonsf = w[fI]*sf[fI];
        const vector wneisf = (1-w[fI])*sf[fI];
        if (markCell[ownFI])
        {
            u[fI] = invPseudoVol[ownFI] & wneisf*vol[ownFI];
            d[ownFI] += invPseudoVol[ownFI] & wonsf*vol[ownFI];
        }
        else
        {
            u[fI] = wneisf;
            d[ownFI] += wonsf;

            if (deferred)
            {
                source[ownFI] +=
                    (wneisf-(invPseudoVol[ownFI] & wneisf*vol[ownFI]))
                   *vf[neiFI];
                source[ownFI] +=
                    (wonsf-(invPseudoVol[ownFI] & wonsf*vol[ownFI]))
                   *vf[ownFI];
            }
        }
        if (markCell[neiFI])
        {
            l[fI] = -invPseudoVol[neiFI] & wonsf*vol[neiFI];
            d[neiFI] -= invPseudoVol[neiFI] & wneisf*vol[neiFI];
        }
        else
        {
            l[fI] = -wonsf;
            d[neiFI] -= wneisf;

            if (deferred)
            {
                source[neiFI] -=
                    (wonsf-(invPseudoVol[neiFI] & wonsf*vol[neiFI]))
                   *vf[ownFI];
                source[neiFI] -=
                    (wonsf-(invPseudoVol[neiFI] & wneisf*vol[neiFI]))
                   *vf[neiFI];
            }
        }
    }

    // Boundary contributions
    forAll(vf.boundaryField(), pI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = patch.Sf();
        const fvsPatchScalarField& pw = tweights().boundaryField()[pI];
        const labelList& fc = patch.faceCells();

        const scalarField iCoeffs(pf.valueInternalCoeffs(pw));
        const scalarField bCoeffs(pf.valueBoundaryCoeffs(pw));

        // Diag contribution
        forAll(pf, fI)
        {
            const label fcI = fc[fI];
            const vector icoeffSfI = iCoeffs[fI]*pSf[fI];
            if (markCell[fcI])
            {
                d[fcI] += invPseudoVol[fcI]&icoeffSfI*vol[fcI];
            }
            else
            {
                d[fcI] += icoeffSfI;
            }
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                bs.coupleUpper()[pI].asLinear();
            CoeffField<vector>::linearTypeField& pcLower =
                bs.coupleLower()[pI].asLinear();

            const vectorField pcl (iCoeffs*pSf);
            const vectorField pcu (-bCoeffs*pSf);

            forAll(pcUpper, fI)
            {
                const label fcI = fc[fI];
                if (markCell[fcI])
                {
                    const scalar& volI = vol[fcI];
                    pcUpper[fI] += invPseudoVol[fcI]&pcu[fI]*volI;
                    pcLower[fI] += invPseudoVol[fcI]&pcl[fI]*volI;
                }
                else
                {
                    pcUpper[fI] += pcu[fI];
                    pcLower[fI] += pcl[fI];
                }
            }
        }
        else
        {
            forAll(pf, fI)
            {
                const label fcI = fc[fI];
                const vector bcoeffSfI = bCoeffs[fI]*pSf[fI];
                if (markCell[fcI])
                {
                    source[fcI] -= invPseudoVol[fcI]&
                        bcoeffSfI*vol[fcI];
                }
                else
                {
                    source[fcI] -= bcoeffSfI;
                    if (deferred)
                    {
                        const vector icoeffSfI = iCoeffs[fI]*pSf[fI];
                        source[fcI] +=
                            (icoeffSfI-(invPseudoVol[fcI]&icoeffSfI*vol[fcI]))
                           *vf[fcI];
                        source[fcI] +=
                            (bcoeffSfI-(invPseudoVol[fcI]&bcoeffSfI*vol[fcI]));
                    }
                }
            }
        }
    }

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

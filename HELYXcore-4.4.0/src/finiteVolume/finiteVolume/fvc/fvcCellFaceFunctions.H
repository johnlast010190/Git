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
    (c) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    (c) 2014-2018 Johan Heyns - CSIR, South Africa
    (c) 2021-2024 Engys Ltd.

Description
    Handy low-level functions to avoid writing face loops

\*---------------------------------------------------------------------------*/

#ifndef fvcCellFaceFunctions_H
#define fvcCellFaceFunctions_H

#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{

template<class Type>
tmp<SurfaceField<Type>>
ownerNeighbourSum
(
    const VolField<Type>& vf,
    const tmp<surfaceScalarField>& tsfP,
    const tmp<surfaceScalarField>& tsfN
)
{
    const surfaceScalarField& sfP = tsfP();
    const surfaceScalarField& sfN = tsfN();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& sfPi = sfP.internalField();
    const scalarField& sfNi = sfN.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            vf.name() + '*' + sfP.name() + '+' + vf.name() + '*' + sfN.name(),
            vf.db(),
            mesh,
            vf.dimensions()*sfP.dimensions()
        )
    );
    SurfaceField<Type>& sf = tsf();

    if (sfP.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction
            << nl << "Dimensions of sfP and sfN differ: " << nl
            << sfP.dimensions() << " != " << sfN.dimensions() << nl
            << exit(FatalError);
    }

    Field<Type>& sfi = sf.internalField();

    forAll(P, fi)
    {
        sfi[fi] = sfPi[fi]*vfi[P[fi]] + sfNi[fi]*vfi[N[fi]];
    }

    forAll(sfP.boundaryField(), pi)
    {
        const fvsPatchScalarField& psfP = sfP.boundaryField()[pi];
        const fvsPatchScalarField& psfN = sfN.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryField()[pi] =
                psfP*vf.boundaryField()[pi].patchInternalField()
              + psfN*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = (psfP+psfN)*vf.boundaryField()[pi];
        }
    }

    tsfP.clear();
    tsfN.clear();

    return tsf;
}


template<class Type>
tmp<Field<Type>> ownerField(const Field<Type>& vfi, const polyMesh& mesh)
{
    const labelUList& P = mesh.faceOwner();
    tmp<Field<Type>> tsf(new Field<Type>(P.size()));
    Field<Type>& sf = tsf.ref();

    forAll(P, fi)
    {
        sf[fi] = vfi[P[fi]];
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> ownerField(const VolField<Type>& vf)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            "owner(" + vf.name() + ')',
            vf.db(),
            mesh,
            vf.dimensions()
        )
    );
    SurfaceField<Type>& sf = tsf.ref();

    forAll(P, fi)
    {
        sf[fi] = vfi[P[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
        sf.boundaryFieldRef()[pi] =
            vf.boundaryField()[pi].patchInternalField();
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> ownerField(const tmp<VolField<Type>> tvf)
{
    tmp<SurfaceField<Type>> tsf = ownerField(tvf());
    tvf.clear();
    return tsf;
}


template<class Type>
tmp<Field<Type>> neighbourField(const Field<Type>& vfi, const polyMesh& mesh)
{
    const labelUList& N = mesh.faceNeighbour();
    tmp<Field<Type>> tsf(new Field<Type>(N.size()));
    Field<Type>& sf = tsf.ref();

    forAll(N, fi)
    {
        sf[fi] = vfi[N[fi]];
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> neighbourField(const VolField<Type>& vf)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& N = mesh.neighbour();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            "neighbour(" + vf.name() + ')',
            vf.db(),
            mesh,
            vf.dimensions()
        )
    );
    SurfaceField<Type>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    forAll(N, fi)
    {
        sfi[fi] = vfi[N[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryFieldRef()[pi] =
                vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> neighbourField(const tmp<VolField<Type>> tvf)
{
    tmp<SurfaceField<Type>> tsf = neighbourField(tvf());
    tvf.clear();
    return tsf;
}


template<class Type>
tmp<VolField<Type>> surfaceSum
(
    const SurfaceField<Type>& sfO,
    const SurfaceField<Type>& sfN
)
{
    const fvMesh& mesh = sfO.mesh();

    tmp<VolField<Type>> tvf
    (
        VolField<Type>::New
        (
            "surfaceSum(" + sfO.name() + ","+  sfN.name()+')',
            sfO.db(),
            mesh,
            dimensioned<Type>("0", sfO.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    VolField<Type>& vf = tvf();

    if (sfO.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction
            << nl << "Dimensions of sfO and sfN differ: " << nl
            << sfO.dimensions() << " != " << sfN.dimensions() << nl
            << exit(FatalError);
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sfO[facei];
        vf[neighbour[facei]] += sfN[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psfO = sfO.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += psfO[facei];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<VolField<Type>>
surfaceSum(tmp<SurfaceField<Type>> tsfO, tmp<SurfaceField<Type>> tsfN)
{
    tmp<VolField<Type>> tvf = surfaceSum(tsfO(), tsfN());
    tsfO.clear();
    tsfN.clear();
    return tvf;
}


template<class Type>
tmp<VolField<Type>>
intSurfaceSum(const SurfaceField<Type>& sfO, const SurfaceField<Type>& sfN)
{
    const fvMesh& mesh = sfO.mesh();

    tmp<VolField<Type>> tvf
    (
        VolField<Type>::New
        (
            "surfaceSum(" + sfO.name() + "," + sfN.name() + ')',
            sfO.db(),
            mesh,
            dimensioned<Type>("0", sfO.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    VolField<Type>& vf = tvf();

    if (sfO.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction
            << nl << "Dimensions of sfO and sfN differ: " << nl
            << sfO.dimensions() << " != " << sfN.dimensions() << nl
            << exit(FatalError);
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sfO[facei];
        vf[neighbour[facei]] += sfN[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psfO = sfO.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                vf[pFaceCells[facei]] += psfO[facei];
            }
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}

template<class Type>
tmp<VolField<Type>>
intSurfaceSum(tmp<SurfaceField<Type>> tsfO, tmp<SurfaceField<Type>> tsfN)
{
    tmp<VolField<Type>> tvf = intSurfaceSum(tsfO(), tsfN());
    tsfO.clear();
    tsfN.clear();
    return tvf;
}


template<class Type>
tmp<VolField<Type>> intSurfaceSum(const SurfaceField<Type>& sf)
{
    const fvMesh& mesh = sf.mesh();

    tmp<VolField<Type>> tvf
    (
        VolField<Type>::New
        (
            "surfaceSum(" + sf.name() + ')',
            sf.db(),
            mesh,
            dimensioned<Type>("0", sf.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    VolField<Type>& vf = tvf();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sf[facei];
        vf[neighbour[facei]] += sf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psf = sf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                vf[pFaceCells[facei]] += psf[facei];
            }
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


// Select from one or the other surface field comparing to a threshold value
template<class Type>
tmp<SurfaceField<Type>>
surfaceFieldSelect
(
    const SurfaceField<Type>& sfL,
    const SurfaceField<Type>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    const fvMesh& mesh = sfL.mesh();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            "select(" + sfL.name() + "," + sfR.name() + ')',
            sfL.db(),
            mesh,
            sfL.dimensions()
        )
    );
    SurfaceField<Type>& sf = tsf.ref();

    if (sfL.dimensions() != sfR.dimensions())
    {
        FatalErrorInFunction
            << "Dimensions of sfL and sfR differ: " << nl
            << sfL.dimensions() << " != " << sfL.dimensions() << nl
            << exit(FatalError);
    }

    forAll(select, facei)
    {
        sf[facei] = (select[facei] >= thres ? sfL[facei] : sfR[facei]);
    }

    forAll(select.boundaryField(), patchi)
    {
        Field<Type>& psf = sf.boundaryFieldRef()[patchi];
        const Field<Type>& psfL = sfL.boundaryField()[patchi];
        const Field<Type>& psfR = sfR.boundaryField()[patchi];
        const scalarField& pSel = select.boundaryField()[patchi];
        forAll(pSel, facei)
        {
            psf[facei] = (pSel[facei] >= thres ? psfL[facei] : psfR[facei]);
        }
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>>
surfaceFieldSelect
(
    tmp<SurfaceField<Type>>& sfL,
    tmp<SurfaceField<Type>>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    tmp<SurfaceField<Type>> tsf =
        surfaceFieldSelect(sfL(), sfR(), select, thres);
    sfL.clear();
    sfR.clear();
    return tsf;
}


template<class Type, class Lambda>
tmp<SurfaceField<Type>>
calcForEachFace(Lambda fn, const SurfaceField<Type>& sf)
{
    const fvMesh& mesh = sf.mesh();

    tmp<SurfaceField<Type>> tsfout
    (
        SurfaceField<Type>::New
        (
            "calc(" + sf.name() + ')',
            sf.db(),
            mesh,
            sf.dimensions()
        )
    );
    SurfaceField<Type>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf = sf.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf[facei]);
        }
    }

    return tsfout;
}


// Call a function for each face of a field; return result as a field
template<class Type, class Lambda>
tmp<SurfaceField<Type>>
calcForEachFace
(
    Lambda fn,
    const SurfaceField<Type>& sf1,
    const SurfaceField<Type>& sf2
)
{
    const fvMesh& mesh = sf1.mesh();

    tmp<SurfaceField<Type>> tsfout
    (
        SurfaceField<Type>::New
        (
            "calc(" + sf1.name() + ',' + sf2.name() + ')',
            sf1.db(),
            mesh,
            sf1.dimensions()
        )
    );
    SurfaceField<Type>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf1[facei], sf2[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf1 = sf1.boundaryField()[patchi];
        const Field<Type>& psf2 = sf2.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf1[facei], psf2[facei]);
        }
    }

    return tsfout;
}


// Sum the co-incident internal and GIB face contributions and sync both to the
// owner and neighbour fields
template<class Type>
void sumIndirectAndInternalFaces
(
    SurfaceField<Type>& sfOwn,
    SurfaceField<Type>& sfNei
)
{
    typename SurfaceField<Type>::Boundary& bsfOwn = sfOwn.boundaryFieldRef();
    typename SurfaceField<Type>::Boundary& bsfNei = sfNei.boundaryFieldRef();

    const fvMesh& mesh = sfOwn.mesh();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            const indirectPolyPatch& ipp =
                refCast<const indirectPolyPatch>(fvP.patch());
            const labelList& addr = ipp.addressing();
            const boolList& fm = ipp.fm();

            fvsPatchField<Type>& psfOwn = bsfOwn[patchi];
            fvsPatchField<Type>& psfNei = bsfNei[patchi];

            if (fvP.isActive())
            {
                bool master = (ipp.indirectPolyPatchType() == "master");
                forAll(addr, pfacei)
                {
                    const label& fi = addr[pfacei];
                    if (mesh.isInternalFace(fi))
                    {
                        if (master != fm[pfacei])
                        {
                            Type sum = sfOwn[fi] + psfOwn[pfacei];
                            sfOwn[fi] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                        else
                        {
                            Type sum = sfNei[fi] + psfOwn[pfacei];
                            sfNei[fi] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                    }
                    else
                    {
                        label pI = mesh.boundaryMesh().whichPatch(fi);
                        label lpfI = mesh.boundaryMesh()[pI].whichFace(fi);
                        if (master != fm[fi])
                        {
                            Type sum = bsfOwn[pI][lpfI] + psfOwn[pfacei];
                            bsfOwn[pI][lpfI] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                        if (master != fm[fi])
                        {
                            Type sum = bsfNei[pI][lpfI] + psfOwn[pfacei];
                            bsfNei[pI][lpfI] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void applyFaceMaskTo(SurfaceField<Type>& sf)
{
    const fvMesh& mesh = sf.mesh();
    Field<Type>& isf = sf.primitiveFieldRef();
    typename SurfaceField<Type>::Boundary& bsf =
        sf.boundaryFieldRef();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            // This blending factor is the proportion of the internal (non indirect-
            // patch) face value that is included for a 'porous' GIB boundary.
            // Currently a placeholder - still to be made settable and generalised
            // to a per-patch/per-face basis
            scalar blendingFactor = 0;

            const indirectPolyPatch& ipp =
                refCast<const indirectPolyPatch>(fvP.patch());
            const labelList& addr = ipp.addressing();

            if (fvP.isActive())
            {
                bool master = (ipp.indirectPolyPatchType() == "master");
                forAll(addr, pfacei)
                {
                    const label& fi = addr[pfacei];
                    if (mesh.isInternalFace(fi))
                    {
                        if (master)
                        {
                            isf[fi] *= blendingFactor;
                        }
                    }
                    else
                    {
                        label pI = mesh.boundaryMesh().whichPatch(fi);
                        label lpfI = mesh.boundaryMesh()[pI].whichFace(fi);
                        bsf[pI][lpfI] *= blendingFactor;
                    }
                }

                bsf[patchi] *= scalar(1)-blendingFactor;
            }
            else
            {
                bsf[patchi] *= Zero;
            }
        }
    }
}


template<class Type>
tmp<SurfaceField<Type>> applyFaceMask(const SurfaceField<Type>& sf)
{
    const fvMesh& mesh = sf.mesh();
    bool haveIndirectPatches = false;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            haveIndirectPatches = true;
            break;
        }
    }

    if (!haveIndirectPatches)
    {
        // Return a const ref; no copy
        return tmp<SurfaceField<Type>>(sf);
    }

    typedef
        SlicedGeometricField
        <
            Type,
            fvsPatchField,
            slicedFvsPatchField,
            surfaceMesh
        >
        slicedFieldType;
    if (isA<slicedFieldType>(sf))
    {
        // For a sliced field, we need to force a deep copy, otherwise the
        // original referenced field will be altered
        tmp<SurfaceField<Type>> tmsf
        (
            new SurfaceField<Type>(sf, sf.mesh(), sf.dimensions())
        );
        tmsf.ref() = sf;
        applyFaceMaskTo(tmsf.ref());
        return tmsf;
    }
    else
    {
        tmp<SurfaceField<Type>> tmsf(new SurfaceField<Type>(sf));
        applyFaceMaskTo(tmsf.ref());
        return tmsf;
    }
}


template<class Type>
tmp<SurfaceField<Type>> applyFaceMask(const tmp<SurfaceField<Type>>& tsf)
{
    applyFaceMaskTo(tsf.ref());
    return tsf;
}


}  // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //


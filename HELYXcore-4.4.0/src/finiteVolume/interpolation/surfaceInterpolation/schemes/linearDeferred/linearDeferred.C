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
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/schemes/linearDeferred/linearDeferred.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::linearDeferred<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "linearDeferred::correction(" + vf.name() + ')',
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    SurfaceField<Type>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const surfaceScalarField& cdWeights
    (
        this->mesh().surfaceInterpolation::weights()
    );

    //correction is linear minus upwind
    forAll(faceFlux, facei)
    {
        scalar fw(cdWeights[facei]);
        // limit weights for stability
        fw = max(0.1, min(fw, 0.9));
        label celli = (faceFlux[facei] >= 0) ? owner[facei] : neighbour[facei];
        sfCorr[facei]
            = fw*vf[owner[facei]] + (1-fw)*vf[neighbour[facei]] - vf[celli];
    }


    typename SurfaceField<Type>::Boundary&
        bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner =
                mesh.boundary()[patchi].faceCells();

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const Field<Type> vfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()()
            );

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                // limit weights for stability
                scalar fw(cdWeights.boundaryField()[patchi][facei]);
                fw = max(0.1, min(fw, 0.9));

                pSfCorr[facei] = fw*vf[own] + (1-fw)*vfNei[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] -= vf[own];
                }
                else
                {
                    pSfCorr[facei] -= vfNei[facei];
                }
            }
        }
    }

    return tsfCorr;
}


namespace Foam
{
    //makelimitedSurfaceInterpolationScheme(linearDeferred)
    makelimitedSurfaceInterpolationTypeScheme(linearDeferred, scalar)
    makelimitedSurfaceInterpolationTypeScheme(linearDeferred, vector)
}

// ************************************************************************* //

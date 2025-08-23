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
    (c) 2010-2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/schemes/linearUpwind/linearUpwind.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::linearUpwind<Type>::correction(const VolField<Type>& vf) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "linearUpwind::correction(" + vf.name() + ')',
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero),
            fvsPatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );

    SurfaceField<Type>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<fv::gradScheme<scalar>> gradScheme_
    (
        fv::gradScheme<scalar>::New
        (
            mesh,
            vf.db(),
            mesh.schemes().gradScheme(gradSchemeName_)
        )
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tmp<volVectorField> tgradVf =
            gradScheme_().grad(vf.component(cmpt), gradSchemeName_);

        const volVectorField& gradVf = tgradVf();

        forAll(faceFlux, facei)
        {
            const label celli =
                (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];

            setComponent(sfCorr[facei], cmpt) =
                (Cf[facei] - C[celli]) & gradVf[celli];
        }

        typename SurfaceField<Type>::
            Boundary& bSfCorr = sfCorr.boundaryFieldRef();

        forAll(bSfCorr, patchi)
        {
            fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

            if (pSfCorr.coupled())
            {
                const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
                const vectorField& pCf = Cf.boundaryField()[patchi];
                const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

                const vectorField pGradVfNei
                (
                    gradVf.boundaryField()[patchi].patchNeighbourField()
                );

                // Build the d-vectors
                const vectorField pd
                (
                    Cf.boundaryField()[patchi].patch().delta()
                );

                forAll(pOwner, facei)
                {
                    label own = pOwner[facei];

                    if (pFaceFlux[facei] > 0)
                    {
                        setComponent(pSfCorr[facei], cmpt) =
                            (pCf[facei] - C[own])
                          & gradVf[own];
                    }
                    else
                    {
                        setComponent(pSfCorr[facei], cmpt) =
                            (pCf[facei] - pd[facei] - C[own])
                          & pGradVfNei[facei];
                    }
                }
            }
        }
    }

    return tsfCorr;
}


template<>
Foam::tmp<Foam::surfaceVectorField>
Foam::linearUpwind<Foam::vector>::correction
(
    const volVectorField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceVectorField> tsfCorr
    (
        surfaceVectorField::New
        (
            "linearUpwind::correction(" + vf.name() + ')',
            mesh,
            dimensioned<vector>(vf.name(), vf.dimensions(), Zero)
        )
    );

    surfaceVectorField& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<fv::gradScheme<vector>> gradScheme_
    (
        fv::gradScheme<vector>::New
        (
            mesh,
            vf.db(),
            mesh.schemes().gradScheme(gradSchemeName_)
        )
    );

    tmp<volTensorField> tgradVf = gradScheme_().grad(vf, gradSchemeName_);
    const volTensorField& gradVf = tgradVf();

    forAll(faceFlux, facei)
    {
        const label celli =
            (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];
        sfCorr[facei] = (Cf[facei] - C[celli]) & gradVf[celli];
    }


    typename surfaceVectorField::Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchVectorField& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = Cf.boundaryField()[patchi];
            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const tensorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
                }
                else
                {
                    pSfCorr[facei] =
                        (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];
                }
            }
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makelimitedSurfaceInterpolationScheme(linearUpwind)
}

// ************************************************************************* //

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
    (c) 2017 Wikki Ltd
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/limitedSchemes/reconCentral/reconCentral.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::reconCentral<Type>::interpolate(const VolField<Type>& vf) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            "interpolate("+vf.name()+')',
            vf.db(),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    SurfaceField<Type>& sf = tsf.ref();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    VolField<typename outerProduct<vector, Type>::type> gradVf
    (
        gradScheme_().grad(vf)
    );

    // Since grad is used on coupled boundaries, correctBoundaryConditions
    // needs to be called.  HJ, 1/Nov/2012
    gradVf.correctBoundaryConditions();

    Field<Type>& sfIn = sf.primitiveFieldRef();

    forAll(sfIn, facei)
    {
        // Owner contribution
        label own = owner[facei];
        sfIn[facei] += 0.5*(vf[own] + ((Cf[facei] - C[own]) & gradVf[own]));

        // Neighbour contribution
        label nei = neighbour[facei];
        sfIn[facei] += 0.5*(vf[nei] + ((Cf[facei] - C[nei]) & gradVf[nei]));
    }


    typename SurfaceField<Type>::Boundary& bSf = sf.boundaryFieldRef();

    forAll(bSf, patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];

        fvsPatchField<Type>& pSf = bSf[patchi];

        const labelUList& pOwner = p.faceCells();

        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (pSf.coupled())
        {
            Field<Type> vfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()
            );

            Field<typename outerProduct<vector, Type>::type> pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors.  Used to calculate neighbour face centre
            // HJ, 19/Apr/2010
            // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
            vectorField pd( p.delta() );

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                // Owner contribution
                pSf[facei] +=
                    0.5*(vf[own] + ((pCf[facei] - C[own]) & gradVf[own]));

                // Neighbour contribution
                pSf[facei] +=
                    0.5*
                    (
                        vfNei[facei]
                      + (
                            (pCf[facei] - pd[facei] - C[own])
                          & pGradVfNei[facei]
                        )
                    );
            }
        }
        else if (vf.boundaryField()[patchi].fixesValue())
        {
            // For fixed boundary patches copy the value
            pSf = vf.boundaryField()[patchi];
        }
        else
        {
            // For patches that do not fix the value, calculate
            // extrapolated field
            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                pSf[facei] =
                    (vf[own] + ((pCf[facei] - C[own]) & gradVf[own]));
            }
        }
    }

    return tsf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::reconCentral<Type>::correction(const VolField<Type>& vf) const
{
    // Note: Correction is calculated by assembling the complete interpolation
    // including extrapolated gradient contribution and subtracting the
    // implicit contribution.  HJ, 27/Mar/2010
    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "reconCentralCorrection(" + vf.name() + ')',
            reconCentral<Type>::interpolate(vf)
          - surfaceInterpolationScheme<Type>::interpolate
            (
                vf,
                this->weights()
            )
        )
    );

    return tsfCorr;
}


namespace Foam
{
    //makelimitedSurfaceInterpolationScheme(reconCentral)
    makelimitedSurfaceInterpolationTypeScheme(reconCentral, scalar)
    makelimitedSurfaceInterpolationTypeScheme(reconCentral, vector)
}

// ************************************************************************* //

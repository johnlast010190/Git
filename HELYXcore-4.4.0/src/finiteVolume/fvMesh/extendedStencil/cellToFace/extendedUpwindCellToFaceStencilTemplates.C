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

\*---------------------------------------------------------------------------*/

#include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::extendedUpwindCellToFaceStencil::weightedSum
(
    const surfaceScalarField& phi,
    const VolField<Type>& fld,
    const List<List<scalar>>& ownWeights,
    const List<List<scalar>>& neiWeights
) const
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> ownFld;
    collectData(ownMap(), ownStencil(), fld, ownFld);
    List<List<Type>> neiFld;
    collectData(neiMap(), neiStencil(), fld, neiFld);

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            fld.name(),
            mesh,
            mesh,
            dimensioned<Type>(fld.name(), fld.dimensions(), Zero)
        )
    );
    SurfaceField<Type>& sf = tsfCorr.ref();

    // Internal faces
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (phi[facei] > 0)
        {
            // Flux out of owner. Use upwind (= owner side) stencil.
            const List<Type>& stField = ownFld[facei];
            const List<scalar>& stWeight = ownWeights[facei];

            forAll(stField, i)
            {
                sf[facei] += stField[i]*stWeight[i];
            }
        }
        else
        {
            const List<Type>& stField = neiFld[facei];
            const List<scalar>& stWeight = neiWeights[facei];

            forAll(stField, i)
            {
                sf[facei] += stField[i]*stWeight[i];
            }
        }
    }

    // Boundaries. Either constrained or calculated so assign value
    // directly (instead of nicely using operator==)
    typename SurfaceField<Type>::
        Boundary& bSfCorr = sf.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            label facei = pSfCorr.patch().start();

            forAll(pSfCorr, i)
            {
                if (phi.boundaryField()[patchi][i] > 0)
                {
                    // Flux out of owner. Use upwind (= owner side) stencil.
                    const List<Type>& stField = ownFld[facei];
                    const List<scalar>& stWeight = ownWeights[facei];

                    forAll(stField, j)
                    {
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                }
                else
                {
                    const List<Type>& stField = neiFld[facei];
                    const List<scalar>& stWeight = neiWeights[facei];

                    forAll(stField, j)
                    {
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                }
                facei++;
            }
        }
    }

    return tsfCorr;
}


// ************************************************************************* //

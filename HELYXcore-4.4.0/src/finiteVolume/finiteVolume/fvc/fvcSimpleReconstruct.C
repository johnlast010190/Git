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
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcReconstruct.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<typename outerProduct<vector,Type>::type>> reconstruct
(
    const SurfaceField<Type>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<SurfaceField<Type>> tmssf = fvc::applyFaceMask(ssf);
    const SurfaceField<Type>& mssf = tmssf();

    tmp<VolField<GradType>> treconField
    (
        VolField<GradType>::New
        (
            "reconstruct(" + ssf.name() + ')',
            ssf.db(),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimArea,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );

    Field<GradType>& rf = treconField();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        rf[own] += (Cf[facei] - C[own])*mssf[facei];
        rf[nei] -= (Cf[facei] - C[nei])*mssf[facei];
    }

    const typename SurfaceField<Type>::
    Boundary& bmsf = mssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchField<Type>& pmsf = bmsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];
            rf[own] += (pCf[pFacei] - C[own])*pmsf[pFacei];
        }
    }

    rf /= mesh.V();

    treconField().correctBoundaryConditions();

    return treconField;
}


template<class Type>
tmp<VolField<typename outerProduct<vector, Type>::type>> reconstruct
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<VolField<GradType>> tvf(fvc::reconstruct(tssf()));
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

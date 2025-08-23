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
    (c) 2017-2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
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
void surfaceIntegrate(Field<Type>& ivf, const SurfaceField<Type>& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    tmp<SurfaceField<Type>> tmssf =
        fvc::applyFaceMask(ssf);
    const Field<Type>& imssf = tmssf();

    forAll(owner, facei)
    {
        ivf[owner[facei]] += imssf[facei];
        ivf[neighbour[facei]] -= imssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pmssf = tmssf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pmssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}


template<class Type>
tmp<VolField<Type>> surfaceIntegrate(const SurfaceField<Type>& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> tvf
    (
        VolField<Type>::New
        (
            "surfaceIntegrate(" + ssf.name() + ')',
            ssf.db(),
            mesh,
            dimensioned<Type>("0", ssf.dimensions()/dimVol, Zero),
            extrapolatedCalculatedFvPatchField<Type>::typeName,
            ssf.boundaryField().patchTypes()
        )
    );
    VolField<Type>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf);

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<VolField<Type>> surfaceIntegrate(const tmp<SurfaceField<Type>>& tssf)
{
    tmp<VolField<Type>> tvf(fvc::surfaceIntegrate(tssf()));
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<VolField<Type>>surfaceSum(const SurfaceField<Type>& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> tvf
    (
        VolField<Type>::New
        (
            "surfaceSum(" + ssf.name() + ')',
            ssf.db(),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), Zero),
            extrapolatedCalculatedFvPatchField<Type>::typeName,
            ssf.boundaryField().patchTypes()
        )
    );
    VolField<Type>& vf = tvf.ref();

    tmp<SurfaceField<Type>> tmssf = fvc::applyFaceMask(ssf);
    const Field<Type>& imssf = tmssf();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += imssf[facei];
        vf[neighbour[facei]] += imssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pmssf = tmssf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += pmssf[facei];
        }
    }

    vf.correctBoundaryConditions();

    //check for zero weights on coupled boundaries and set to zero grad value
    //addresses issue with un-covered AMI faces and fvc::reconstruct
    typename VolField<Type>::Boundary& vfbf = vf.boundaryFieldRef();
    forAll(vfbf, patchi)
    {
        if (vfbf[patchi].coupled())
        {
            Field<Type> pif(vfbf[patchi].patchInternalField());

            forAll(vfbf[patchi], fi)
            {
                if (vfbf[patchi].patch().weights()[fi] == 0)
                {
                    vfbf[patchi][fi] = pif[fi];
                }
            }
        }
    }

    return tvf;
}


template<class Type>
tmp<VolField<Type>> surfaceSum(const tmp<SurfaceField<Type>>& tssf)
{
    tmp<VolField<Type>> tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

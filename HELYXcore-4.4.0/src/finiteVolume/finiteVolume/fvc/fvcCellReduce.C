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

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcCellReduce.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
Foam::tmp<Foam::VolField<Type>> Foam::fvc::cellReduce
(
    const SurfaceField<Type>& ssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> tresult
    (
        VolField<Type>::New
        (
            "cellReduce(" + ssf.name() + ')',
            ssf.db(),
            mesh,
            dimensioned<Type>("initialValue", ssf.dimensions(), nullValue),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );

    VolField<Type>& result = tresult.ref();

    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    forAll(own, i)
    {
        const label celli = own[i];
        cop(result[celli], ssf[i]);
    }
    forAll(nbr, i)
    {
        const label celli = nbr[i];
        cop(result[celli], ssf[i]);
    }

    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolField<Type>> Foam::fvc::cellReduce
(
    const tmp<SurfaceField<Type>>& tssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    tmp<VolField<Type>> tvf(cellReduce(tssf, cop, nullValue));
    tssf.clear();

    return tvf;
}


// ************************************************************************* //

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

#include "finiteVolume/fvc/fvcAverage.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fvMesh/fvMesh.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> average(const SurfaceField<Type>& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> taverage
    (
        VolField<Type>::New
        (
            "average(" + ssf.name() + ')',
            ssf.db(),
            mesh,
            ssf.dimensions()
        )
    );
    VolField<Type>& av = taverage.ref();

    av.primitiveFieldRef() =
    (
        surfaceSum(mesh.magSf()*ssf)().primitiveField()
       /surfaceSum(mesh.magSf())().primitiveField()
    );

    typename VolField<Type>::
    Boundary& bav = av.boundaryFieldRef();

    forAll(bav, patchi)
    {
        bav[patchi] = ssf.boundaryField()[patchi];
    }

    av.correctBoundaryConditions();

    return taverage;
}


template<class Type>
tmp<VolField<Type>> average(const tmp<SurfaceField<Type>>& tssf)
{
    tmp<VolField<Type>> taverage(fvc::average(tssf()));
    tssf.clear();
    return taverage;
}


template<class Type>
tmp<VolField<Type>> average(const VolField<Type>& vtf)
{
    return fvc::average(linearInterpolate(vtf));
}


template<class Type>
tmp<VolField<Type>> average(const tmp<VolField<Type>>& tvtf)
{
    tmp<VolField<Type>> taverage(fvc::average(tvtf()));
    tvtf.clear();
    return taverage;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

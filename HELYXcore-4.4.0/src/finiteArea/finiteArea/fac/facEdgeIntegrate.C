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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Engys Ltd.

Description


\*---------------------------------------------------------------------------*/

#include "finiteArea/fac/facEdgeIntegrate.H"
#include "faMesh/faMesh.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
edgeIntegrate
(
    const GeometricField<Type, faePatchField, faEdgeMesh>& ssf
)
{
    const faMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh>> tvf
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "edgeIntegrate("+ssf.name()+')',
                ssf.instance(),
                ssf.db()
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimArea,
                pTraits<Type>::zero
            ),
            zeroGradientFaPatchField<Type>::typeName
        )
    );
    GeometricField<Type, faPatchField, areaMesh>& vf = tvf.ref();


    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] -= ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const labelUList& pEdgeFaces =
            mesh.boundary()[patchI].edgeFaces();

        const faePatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pEdgeFaces[faceI]] += pssf[faceI];
        }
    }

    vf.primitiveFieldRef() /= mesh.S();
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
edgeIntegrate
(
    const tmp<GeometricField<Type, faePatchField, faEdgeMesh>>& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tvf
    (
        fac::edgeIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
edgeSum
(
    const GeometricField<Type, faePatchField, faEdgeMesh>& ssf
)
{
    const faMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh>> tvf
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "edgeSum("+ssf.name()+')',
                ssf.instance(),
                ssf.db()
            ),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
            zeroGradientFaPatchField<Type>::typeName
        )
    );
    GeometricField<Type, faPatchField, areaMesh>& vf = tvf.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] += ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const labelUList& pEdgeFaces =
            mesh.boundary()[patchI].edgeFaces();

        const faePatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pEdgeFaces[faceI]] += pssf[faceI];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>> edgeSum
(
    const tmp<GeometricField<Type, faePatchField, faEdgeMesh>>& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tvf =
        edgeSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

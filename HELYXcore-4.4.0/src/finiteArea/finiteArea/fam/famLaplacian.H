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

Namespace
    famLaplacian

Description
    Calculate the matrix for the laplacian of the field.

SourceFiles
    famLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef famLaplacian_H
#define famLaplacian_H

#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "faMatrices/faMatrix/faMatrix.H"
#include "primitives/zero/zero.H"
#include "primitives/one/one.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fam functions Declaration
\*---------------------------------------------------------------------------*/

namespace fam
{
    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const zero&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const zero&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const one&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const one&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<GType, faPatchField, areaMesh>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<GType, faPatchField, areaMesh>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<GeometricField<GType, faPatchField, areaMesh>>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<GeometricField<GType, faPatchField, areaMesh>>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<GType, faePatchField, faEdgeMesh>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<GType, faePatchField, faEdgeMesh>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type, class GType>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );
    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/fam/famLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

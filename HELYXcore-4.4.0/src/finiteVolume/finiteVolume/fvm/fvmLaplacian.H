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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019 Engys Ltd.

InNamespace
    Foam::fvm

Description
    Calculate the matrix for the laplacian of the field.

SourceFiles
    fvmLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef fvmLaplacian_H
#define fvmLaplacian_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/zero/zero.H"
#include "primitives/one/one.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvMatrix;

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const zero&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const zero&,
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const one&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp
    <
        fvBlockMatrix<Type>
    > bLaplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvm/fvmLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

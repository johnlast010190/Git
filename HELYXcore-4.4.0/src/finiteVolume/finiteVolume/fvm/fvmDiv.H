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
    (c) 2019-2021 Engys Ltd.

InNamespace
    Foam::fvm

Description
    Calculate the matrix for the divergence of the given field and flux.

SourceFiles
    fvmDiv.C

\*---------------------------------------------------------------------------*/

#ifndef fvmDiv_H
#define fvmDiv_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fvMatrices/fvMatrices.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const surfaceScalarField&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&,
        const word& name
    );


    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const surfaceScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const surfaceScalarField&,
        const volScalarField&,
        const VolField<Type>&,
        const word& name=word::null
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const volScalarField&,
        const VolField<Type>&,
        const word& name=word::null
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const surfaceScalarField&,
        const tmp<volScalarField>&,
        const VolField<Type>&,
        const word& name=word::null
    );

    template<class Type>
    tmp<fvMatrix<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const tmp<volScalarField>&,
        const VolField<Type>&,
        const word& name=word::null
    );

    template<class Type>
    tmp<fvBlockMatrix<Type>> bDiv
    (
        const surfaceScalarField&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<fvBlockMatrix<Type>> bDiv
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&,
        const word& name
    );


    template<class Type>
    tmp<fvBlockMatrix<Type>> bDiv
    (
        const surfaceScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<fvBlockMatrix<Type>> bDiv
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > div
    (
        const VolField<Type>&
    );


    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > div
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > div
    (
        const VolField<scalar>&,
        const VolField<Type>&
    );


    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > div
    (
        const VolField<scalar>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > bcDiv
    (
        const SurfaceField<scalar>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > bcDiv
    (
        const tmp<SurfaceField<scalar>>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    > bcDiv
    (
        const SurfaceField<scalar>&,
        const VolField<Type>&,
        const word&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvm/fvmDiv.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

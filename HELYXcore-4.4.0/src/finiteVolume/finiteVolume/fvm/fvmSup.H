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

InNamespace
    Foam::fvm

Description
    Calculate the matrix for implicit and explicit sources.

SourceFiles
    fvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef fvmSup_H
#define fvmSup_H

#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "fields/Fields/zeroField/zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvMatrix;

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
   // Explicit source

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const DimensionedField<Type, volMesh>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<DimensionedField<Type, volMesh>>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<VolField<Type>>&,
            const VolField<Type>&
        );

        template<class Type>
        zeroField Su
        (
            const zero&,
            const VolField<Type>&
        );


    // Implicit source

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const volScalarField::Internal&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<volScalarField::Internal>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );


        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );


        template<class Type>
        zeroField Sp
        (
            const zero&,
            const VolField<Type>&
        );



        template<class Type>
        tmp<fvBlockMatrix<Type>> bSp
        (
            const tmp
            <
                GeometricField
                <
                    typename outerProduct<vector, Type>::type,
                    fvPatchField,
                    volMesh
                >
            >&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvBlockMatrix<Type>> bSuSp
        (
            const tmp
            <
                GeometricField
                <
                    typename outerProduct<vector, Type>::type,
                    fvPatchField,
                    volMesh
                >
            >&,
            const VolField<Type>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const volScalarField::Internal&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<volScalarField::Internal>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );

        template<class Type>
        zeroField SuSp
        (
            const zero&,
            const VolField<Type>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvm/fvmSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

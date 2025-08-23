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
    Foam::fvc

Description
    Calculate the field for explicit evaluation of implicit and explicit
    sources.

SourceFiles
    fvcSup.C

\*---------------------------------------------------------------------------*/

#ifndef fvcSup_H
#define fvcSup_H

#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
   // Explicit source

        template<class Type>
        tmp<VolField<Type>> Su
        (
            const VolField<Type>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolField<Type>> Su
        (
            const tmp<VolField<Type>>&,
            const VolField<Type>&
        );


    // Implicit source

        template<class Type>
        tmp<VolField<Type>> Sp
        (
            const volScalarField&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolField<Type>> Sp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );


        template<class Type>
        tmp<VolField<Type>> Sp
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<VolField<Type>> SuSp
        (
            const volScalarField&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolField<Type>> SuSp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2017 Engys Ltd.

InNamespace
    Foam::fvc

Description
    Reconstruct volField from a face flux field.

SourceFiles
    fvcReconstruct.C

\*---------------------------------------------------------------------------*/


#ifndef fvcReconstruct_H
#define fvcReconstruct_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<typename outerProduct<vector, Type>::type>> reconstruct
    (
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolField<typename outerProduct<vector, Type>::type>> reconstruct
    (
        const tmp<SurfaceField<Type>>&
    );

    tmp<volScalarField> reconstructMag(const surfaceScalarField&);
    tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>&);

    template<class Type>
    tmp<VolField<typename outerProduct<vector, Type>::type>> reconstructV
    (
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolField<typename outerProduct<vector, Type>::type>> reconstructV
    (
        const tmp<SurfaceField<Type>>&
    );


    template<class Type>
    tmp<VolField<Type>> reconstructSmooth
    (
        const VolField<Type>&,
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> reconstructSmooth
    (
        const VolField<Type> &,
        const tmp<SurfaceField<Type>>&
    );

    template<class Type>
    tmp<VolField<Type>> reconstructSmooth(const VolField<Type>&);

    template<class Type>
    tmp<VolField<Type>> reconstructSmooth(const tmp<VolField<Type>>&);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcReconstruct.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

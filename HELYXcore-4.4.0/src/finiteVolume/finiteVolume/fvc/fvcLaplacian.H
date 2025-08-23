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
    Calculate the laplacian of the given field.

SourceFiles
    fvcLaplacian.C

\*---------------------------------------------------------------------------*/


#ifndef fvcLaplacian_H
#define fvcLaplacian_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "dimensionedTypes/dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

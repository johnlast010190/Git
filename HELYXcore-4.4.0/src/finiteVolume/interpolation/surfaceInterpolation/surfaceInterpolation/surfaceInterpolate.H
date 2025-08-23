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
    Surface Interpolation

SourceFiles
    surfaceInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolate_H
#define surfaceInterpolate_H

#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "primitives/one/one.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type>> scheme
    (
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type>> scheme
    (
        const surfaceScalarField& faceFlux,
        const word& name
    );


    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type>> scheme
    (
        const fvMesh& mesh,
        const objectRegistry& db,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type>> scheme
    (
        const fvMesh& mesh,
        const objectRegistry& db,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const tmp<VolField<Type>>& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const tmp<VolField<Type>>& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const tmp<VolField<Type>>& tvf,
        const word& name
    );


    //- Interpolate field onto faces using 'interpolate(\<name\>)'
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& tvf
    );

    //- Interpolate tmp field onto faces using 'interpolate(\<name\>)'
    template<class Type>
    static tmp<SurfaceField<Type>> interpolate
    (
        const tmp<VolField<Type>>& tvf
    );


    //- Interpolate boundary field onto faces (simply a type conversion)
    template<class Type>
    static tmp<FieldField<fvsPatchField, Type>> interpolate
    (
        const FieldField<fvPatchField, Type>& fvpff
    );

    //- Interpolate boundary field onto faces (simply a type conversion)
    template<class Type>
    static tmp<FieldField<fvsPatchField, Type>> interpolate
    (
        const tmp<FieldField<fvPatchField, Type>>& tfvpff
    );

    //- Interpolate 'one' returning 'one'
    inline one interpolate(const one&)
    {
        return one();
    }


    //- Interpolate field onto faces
    //  and 'dot' with given surfaceVectorField Sf
    template<class Type>
    static tmp<SurfaceField<typename innerProduct<vector, Type>::type>>
    dotInterpolate
    (
        const surfaceVectorField& Sf,
        const VolField<Type>& tvf,
        const word = word::null
    );

    //- Interpolate tmp field onto faces
    //  and 'dot' with given surfaceVectorField Sf
    template<class Type>
    static tmp<SurfaceField<typename innerProduct<vector, Type>::type>>
    dotInterpolate
    (
        const surfaceVectorField& Sf,
        const tmp<VolField<Type>>& tvf,
        const word = word::null
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

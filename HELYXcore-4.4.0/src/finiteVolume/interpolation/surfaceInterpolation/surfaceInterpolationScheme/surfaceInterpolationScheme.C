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
    (c) 2010-2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::surfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const objectRegistry& db,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (surfaceInterpolation::debug || surfaceInterpolationScheme<Type>::debug)
    {
        InfoInFunction << "Discretisation scheme = " << schemeName << endl;
    }

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTable_().find(schemeName);

    if (constructorIter == MeshConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, db, schemeData);
}


template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::surfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (surfaceInterpolation::debug || surfaceInterpolationScheme<Type>::debug)
    {
        InfoInFunction
            << "Discretisation scheme = " << schemeName << endl;
    }

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTable_().find(schemeName);

    if (constructorIter == MeshFluxConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::surfaceInterpolationScheme<Type>::~surfaceInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::surfaceInterpolationScheme<Type>::interpolate
(
    const VolField<Type>& vf,
    const tmp<surfaceScalarField>& tlambdas,
    const tmp<surfaceScalarField>& tys
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces "
               "without explicit correction"
            << endl;
    }

    const surfaceScalarField& lambdas = tlambdas();
    const surfaceScalarField& ys = tys();

    const Field<Type>& vfi = vf;
    const scalarField& lambda = lambdas;
    const scalarField& y = ys;

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            "interpolate(" + vf.name() + ')',
            vf.db(),
            mesh,
            vf.dimensions(),
            fvsPatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );
    SurfaceField<Type>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*vfi[P[fi]] + y[fi]*vfi[N[fi]];
    }


    // Interpolate across coupled patches using given lambdas and ys
    typename SurfaceField<Type>::
        Boundary& sfbf = sf.boundaryFieldRef();

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const fvsPatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            sfbf[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
              + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sfbf[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();
    tys.clear();

    return tsf;
}


template<class Type>
template<class SFType>
Foam::tmp
<
    Foam::SurfaceField
    <
        typename Foam::innerProduct<typename SFType::value_type, Type>::type
    >
>
Foam::surfaceInterpolationScheme<Type>::dotInterpolate
(
    const SFType& Sf,
    const VolField<Type>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces "
               "without explicit correction"
            << endl;
    }

    typedef typename Foam::innerProduct<typename SFType::value_type, Type>::type
        RetType;

    const surfaceScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf;
    const scalarField& lambda = lambdas;

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    tmp<SurfaceField<RetType>> tsf
    (
        SurfaceField<RetType>::New
        (
            "interpolate(" + vf.name() + ')',
            vf.db(),
            mesh,
            Sf.dimensions()*vf.dimensions(),
            fvsPatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );
    SurfaceField<RetType>& sf = tsf.ref();

    Field<RetType>& sfi = sf.primitiveFieldRef();

    const typename SFType::Internal& Sfi = Sf();

    for (label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = Sfi[fi] & (lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]]);
    }

    // Interpolate across coupled patches using given lambdas

    typename SurfaceField<RetType>::Boundary& sfbf = sf.boundaryFieldRef();

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const typename SFType::Patch& pSf = Sf.boundaryField()[pi];
        fvsPatchField<RetType>& psf = sfbf[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            psf =
                pSf
              & (
                    pLambda*vf.boundaryField()[pi].patchInternalField()
                  + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField()
                );
        }
        else
        {
            psf = pSf & vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

//    tsf.ref().oriented() = Sf.oriented();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::surfaceInterpolationScheme<Type>::interpolate
(
    const VolField<Type>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    return dotInterpolate(geometricOneField(), vf, tlambdas);
}


template<class Type>
Foam::tmp
<
    Foam::SurfaceField
    <
        typename Foam::innerProduct<Foam::vector, Type>::type
    >
>
Foam::surfaceInterpolationScheme<Type>::dotInterpolate
(
    const surfaceVectorField& Sf,
    const VolField<Type>& vf
) const
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces"
            << endl;
    }

    tmp<SurfaceField<typename Foam::innerProduct<Foam::vector, Type>::type>>
        tsf = dotInterpolate(Sf, vf, weights(vf));

    tsf.ref().oriented() = Sf.oriented();

    if (corrected())
    {
        tsf.ref() += Sf & correction(vf);
    }

    return tsf;
}


template<class Type>
Foam::tmp
<
    Foam::SurfaceField<typename Foam::innerProduct<Foam::vector, Type>::type>
>
Foam::surfaceInterpolationScheme<Type>::dotInterpolate
(
    const surfaceVectorField& Sf,
    const tmp<VolField<Type>>& tvf
) const
{
    tmp<SurfaceField<typename Foam::innerProduct<Foam::vector, Type>::type>>
        tSfDotinterpVf = dotInterpolate(Sf, tvf());

    tvf.clear();
    return tSfDotinterpVf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::surfaceInterpolationScheme<Type>::interpolate
(
    const VolField<Type>& vf
) const
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces"
            << endl;
    }

    tmp<SurfaceField<Type>> tsf = interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::surfaceInterpolationScheme<Type>::interpolate
(
    const tmp<VolField<Type>>& tvf
) const
{
    tmp<SurfaceField<Type>> tinterpVf = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// ************************************************************************* //

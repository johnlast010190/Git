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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.H"
#include "finiteVolume/divSchemes/divScheme/divScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemes().divScheme(name)
    )().fvmDiv(flux, vf);
}

template<class Type>
tmp<fvMatrix<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Div(fvm::div(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvm::div(flux, vf, fvc::divSchemeName(flux, vf));
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> Div(fvm::div(tflux(), vf));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const surfaceScalarField& flux,
    const volScalarField& vfMult,
    const VolField<Type>& vf,
    const word& name
)
{
    if (name == word::null)
    {
        return
            fvm::div
            (
                flux,
                vfMult,
                vf,
                "div("+flux.name()+','+vfMult.name()+'*'+vf.name()+')'
            );
    }
    else
    {
        return fv::convectionScheme<Type>::New
        (
            vf.mesh(),
            flux,
            vf.mesh().schemes().divScheme(name)
        )().fvmDiv(flux, vfMult, vf);
    }
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const volScalarField& vfMult,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Div(fvm::div(tflux(), vfMult, vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const surfaceScalarField& flux,
    const tmp<volScalarField>& tvfMult,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Div(fvm::div(flux, tvfMult(), vf, name));
    tvfMult.clear();
    return Div;
}


template<class Type>
tmp<fvMatrix<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<volScalarField>& tvfMult,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Div(fvm::div(tflux(), tvfMult(), vf, name));
    tflux.clear();
    tvfMult.clear();
    return Div;
}


template<class Type>
tmp<fvBlockMatrix<Type>>
bDiv
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemes().divScheme(name)
    )().fvmBDiv(flux, vf);
}

template<class Type>
tmp<fvBlockMatrix<Type>>
bDiv
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvBlockMatrix<Type>> Div(fvm::bDiv(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<fvBlockMatrix<Type>>
bDiv
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvm::bDiv(flux, vf, fvc::divSchemeName(flux, vf));
}


template<class Type>
tmp<fvBlockMatrix<Type>>
bDiv
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<fvBlockMatrix<Type>> Div(fvm::bDiv(tflux(), vf));
    tflux.clear();
    return Div;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> div
(
    const VolField<Type>& vf
)
{
    return fvm::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> div
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().divScheme(name)
    )().fvmDiv(vf);
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> div
(
    const VolField<scalar>& rho,
    const VolField<Type>& vf
)
{
    return fvm::div(rho, vf, "div("+ rho.name()+','+ vf.name()+')');
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> div
(
    const VolField<scalar>& rho,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().divScheme(name)
    )().fvmDiv(fvc::interpolate(rho)(), vf);
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> bcDiv
(
    const SurfaceField<scalar>& rhof,
    const VolField<Type>& vf
)
{
    return fvm::bcDiv(rhof, vf, "div("+ rhof.name()+','+ vf.name()+')');
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> bcDiv
(
    const tmp<SurfaceField<scalar>>& rhof,
    const VolField<Type>& vf
)
{
    return bcDiv(rhof(), vf);
    rhof.clear();
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> bcDiv
(
    const SurfaceField<scalar>& rhof,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().divScheme(name)
    )().fvmBCDiv(rhof, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

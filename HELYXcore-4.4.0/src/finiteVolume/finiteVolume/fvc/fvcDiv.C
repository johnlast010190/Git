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
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/divSchemes/divScheme/divScheme.H"
#include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> div(const SurfaceField<Type>& ssf)
{
    return VolField<Type>::New
    (
        "div(" + ssf.name() + ')',
        fvc::surfaceIntegrate(ssf)
    );
}


template<class Type>
tmp<VolField<Type>> div(const tmp<SurfaceField<Type>>& tssf)
{
    tmp<VolField<Type>> Div(fvc::div(tssf()));
    tssf.clear();
    return Div;
}


template<class Type>
word divSchemeName
(
    const surfaceScalarField& ssf,
    const VolField<Type>& vf
)
{
    // If the scheme name isn't found, try without the group parts
    const word name1("div(" + ssf.name() + ',' + vf.name() + ')');
    if (!vf.mesh().schemes().divSchemes().found(name1))
    {
        const word name2("div(" + ssf.member() + ',' + vf.member() + ')');
        if (vf.mesh().schemes().divSchemes().found(name2))
        {
            return name2;
        }
    }
    return name1;
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>> div
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(), vf.db(), vf.mesh().schemes().divScheme(name)
    ).ref().fvcDiv(vf);
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>> div
(
    const tmp<VolField<Type>>& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolField<DivType>> Div(fvc::div(tvvf(), name));
    tvvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>> div
(
    const VolField<Type>& vf
)
{
    return fvc::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>> div
(
    const tmp<VolField<Type>>& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolField<DivType>> Div(fvc::div(tvvf()));
    tvvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
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
    ).ref().fvcDiv(flux, vf);
}


template<class Type>
tmp<VolField<Type>> div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolField<Type>> Div(fvc::div(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Div(fvc::div(flux, tvf(), name));
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Div(fvc::div(tflux(), tvf(), name));
    tflux.clear();
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvc::div(flux, vf, divSchemeName(flux, vf));
}


template<class Type>
tmp<VolField<Type>> div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<VolField<Type>> Div(fvc::div(tflux(), vf));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Div(fvc::div(flux, tvf()));
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>> div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Div(fvc::div(tflux(), tvf()));
    tflux.clear();
    tvf.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

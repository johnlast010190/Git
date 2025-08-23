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

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcFlux.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        phi,
        vf.mesh().schemes().divScheme(name)
    )().flux(phi, vf);
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), vf, name)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(phi, tvf(), name)
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), tvf(), name)
    );
    tphi.clear();
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf
)
{
    return fvc::flux
    (
        phi, vf, "flux("+phi.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), vf)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(phi, tvf())
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), tvf())
    );
    tphi.clear();
    tvf.clear();
    return Flux;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

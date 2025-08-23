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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf);
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvc::meshPhi
(
    const VolField<Type>& vf,
    const label patchi
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf, patchi);
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


template<class Type>
void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi -= fvc::meshPhi(vf);
    }
}

template<class Type>
void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi -= rho*fvc::meshPhi(rho, vf);
    }
}

template<class Type>
void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, vf);
    }
}


template<class Type>
void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const surfaceScalarField& rhof,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi -= rhof*fvc::meshPhi(rho, vf);
    }
}


template<class Type>
void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi += fvc::meshPhi(vf);
    }
}

template<class Type>
void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi += rho*fvc::meshPhi(rho, vf);
    }
}

template<class Type>
void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (phi.mesh().moving())
    {
        phi += fvc::interpolate(rho)*fvc::meshPhi(rho, vf);
    }
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fvc::meshPhi(vf);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fvc::interpolate(rho)*fvc::meshPhi(rho, vf);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fvc::meshPhi(vf);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fvc::interpolate(rho)*fvc::meshPhi(rho, vf);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


// ************************************************************************* //

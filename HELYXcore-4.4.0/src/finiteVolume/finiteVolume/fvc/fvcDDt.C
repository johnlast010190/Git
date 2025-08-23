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
    (c) 2020-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcDDt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMesh/fvMesh.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
DDt
(
    const surfaceScalarField& phi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> ddtDivPhiPsi
        = fvc::ddt(psi) + fvc::div(phi, psi);

    if (phi.mesh().moving())
    {
        return ddtDivPhiPsi - fvc::div(fvc::absolute(phi, psi))*psi;
    }
    else
    {
        return ddtDivPhiPsi - fvc::div(phi)*psi;
    }
}

template<class Type>
tmp<VolField<Type>>
DDt
(
    const geometricOneField& rho,
    const surfaceScalarField& phi,
    const VolField<Type>& psi
)
{
    return fvc::DDt(phi, psi);
}

template<class Type>
tmp<VolField<Type>>
DDt
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const VolField<Type>& psi
)
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        return
            1/rho
           *(
               fvc::ddt(rho, psi) + fvc::div(phi, psi)
             - (fvc::ddt(rho) + fvc::div(phi))*psi
            );
    }
    else
    {
        return fvc::DDt(phi, psi);
    }
}

template<class Type>
tmp<VolField<Type>>
DDt
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const geometricOneField& rho,
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(rho, tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const volScalarField& rho,
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(rho, tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const tmp<geometricOneField>& trho,
    const surfaceScalarField& phi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(trho(), phi, psi)
    );
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const tmp<volScalarField>& trho,
    const surfaceScalarField& phi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(trho(), phi, psi)
    );
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const tmp<geometricOneField>& trho,
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(trho(), tphi(), psi)
    );
    tphi.clear();
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<VolField<Type>> DDt
(
    const tmp<volScalarField>& trho,
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& psi
)
{
    tmp<VolField<Type>> DDtPsi
    (
        fvc::DDt(trho(), tphi(), psi)
    );
    tphi.clear();
    trho.clear();
    return DDtPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

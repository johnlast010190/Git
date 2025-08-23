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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> EulerDdtScheme<Type>::fvcDdt(const dimensioned<Type>& dt)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt(" + dt.name() + ')');

    if (mesh().moving())
    {
        tmp<VolField<Type>> tdtdt
        (
            VolField<Type>::New
            (
                ddtName,
                mesh(),
                dimensioned<Type>("0", dt.dimensions()/dimTime, Zero)
            )
        );

        tdtdt.ref().primitiveFieldRef() =
            rDeltaT.value()*dt.value()*(1.0 - mesh().Vsc0()/mesh().Vsc());

        return tdtdt;
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            mesh(),
            dimensioned<Type>("0", dt.dimensions()/dimTime, Zero),
            calculatedFvPatchField<Type>::typeName
        );
    }
}


template<class Type>
tmp<VolField<Type>> EulerDdtScheme<Type>::fvcDdt(const VolField<Type>& vf)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt(" + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT
           *(
                vf()
              - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()
           *(
                vf.boundaryField() - vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New(ddtName, rDeltaT*(vf - vf.oldTime()));
    }
}


template<class Type>
tmp<VolField<Type>> EulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt(" + rho.name() + ',' + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*rho
           *(
                vf()
              - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()*rho.value()
           *(
                vf.boundaryField() - vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New(ddtName, rDeltaT*rho*(vf - vf.oldTime()));
    }
}


template<class Type>
tmp<VolField<Type>> EulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt(" + rho.name() + ',' + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT
           *(
                rho()*vf()
              - rho.oldTime()()
               *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()
           *(
                rho.boundaryField()*vf.boundaryField()
              - rho.oldTime().boundaryField()
               *vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
        );
    }
}


template<class Type>
tmp<VolField<Type>> EulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const word ddtName
    (
        "ddt(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')'
    );

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT
           *(
                alpha()
               *rho()
               *vf()
              - alpha.oldTime()()
               *rho.oldTime()()
               *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()
           *(
                alpha.boundaryField()
               *rho.boundaryField()
               *vf.boundaryField()
              - alpha.oldTime().boundaryField()
               *rho.oldTime().boundaryField()
               *vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT
           *(
               alpha*rho*vf
             - alpha.oldTime()*rho.oldTime()*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<SurfaceField<Type>> EulerDdtScheme<Type>::fvcDdt
(
    const SurfaceField<Type>& sf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    return SurfaceField<Type>::New
    (
        "ddt(" + sf.name() + ')',
        rDeltaT*(sf - sf.oldTime())
    );
}


template<class Type>
tmp<fvMatrix<Type>> EulerDdtScheme<Type>::fvmDdt(const VolField<Type>& vf)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>(vf, vf.dimensions()*dimVol/dimTime)
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> EulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.value()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() =
            rDeltaT*rho.value()*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() =
            rDeltaT*rho.value()*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> EulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() =
            rDeltaT
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() =
            rDeltaT
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> EulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() =
        rDeltaT*alpha.primitiveField()*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() =
            rDeltaT
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() =
            rDeltaT
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtUfCorr
(
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf,
    const word interpolationName
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    const fluxFieldType phiCorr
    (
        phiUf0
      - fvc::dotInterpolate(mesh().Sf(), U.oldTime(), interpolationName)
    );

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr)*rDeltaT*phiCorr
    );
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const VolField<Type>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const fluxFieldType phiCorr
    (
        phi.oldTime()
      - fvc::dotInterpolate(mesh().Sf(), U.oldTime(), interpolationName)
    );

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
       *rDeltaT*phiCorr
    );
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& rhoUf,
    const word interpolationName
)
{
    if
    (
        U.dimensions() == dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

        const VolField<Type> rhoU0(rho.oldTime()*U.oldTime());
        const fluxFieldType phiUf0(mesh().Sf() & rhoUf.oldTime());
        const fluxFieldType phiCorr
        (
            phiUf0
          - fvc::dotInterpolate(mesh().Sf(), rhoU0, interpolationName)
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name()
          + ',' + rhoUf.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr)*rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        return fvcDdtUfCorr(U, rhoUf, interpolationName);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of rhoUf are not correct"
            << abort(FatalError);
        ::abort();
    }
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
        VolField<Type> rhoU0(rho.oldTime()*U.oldTime());

        fluxFieldType phiCorr
        (
            phi.oldTime()
          - fvc::dotInterpolate(mesh().Sf(), rhoU0, interpolationName)
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), phiCorr)
           *rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return fvcDdtPhiCorr(U, phi, interpolationName);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);
        ::abort();
    }
}


template<class Type>
tmp<surfaceScalarField> EulerDdtScheme<Type>::meshPhi(const VolField<Type>&)
{
    return mesh().phi();
}


template<class Type>
tmp<scalarField> EulerDdtScheme<Type>::meshPhi
(
    const VolField<Type>&,
    const label patchi
)
{
    return mesh().phi().boundaryField()[patchi];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

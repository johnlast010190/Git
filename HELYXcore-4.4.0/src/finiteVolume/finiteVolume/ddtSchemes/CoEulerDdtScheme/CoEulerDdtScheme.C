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

#include "finiteVolume/ddtSchemes/CoEulerDdtScheme/CoEulerDdtScheme.H"
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
tmp<volScalarField> CoEulerDdtScheme<Type>::CorDeltaT() const
{
    const surfaceScalarField cofrDeltaT(CofrDeltaT());

    tmp<volScalarField> tcorDeltaT
    (
        volScalarField::New
        (
            "CorDeltaT",
            cofrDeltaT.db(),
            mesh(),
            dimensionedScalar(cofrDeltaT.dimensions(), 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& corDeltaT = tcorDeltaT.ref();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, facei)
    {
        corDeltaT[owner[facei]] =
            max(corDeltaT[owner[facei]], cofrDeltaT[facei]);

        corDeltaT[neighbour[facei]] =
            max(corDeltaT[neighbour[facei]], cofrDeltaT[facei]);
    }

    const surfaceScalarField::Boundary& cofrDeltaTbf =
        cofrDeltaT.boundaryField();

    forAll(cofrDeltaTbf, patchi)
    {
        const fvsPatchScalarField& pcofrDeltaT = cofrDeltaTbf[patchi];
        const fvPatch& p = pcofrDeltaT.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pcofrDeltaT, patchFacei)
        {
            corDeltaT[faceCells[patchFacei]] = max
            (
                corDeltaT[faceCells[patchFacei]],
                pcofrDeltaT[patchFacei]
            );
        }
    }

    corDeltaT.correctBoundaryConditions();

    return tcorDeltaT;
}


template<class Type>
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::CofrDeltaT() const
{
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    const surfaceScalarField& phi =
        static_cast<const objectRegistry&>(mesh())
        .lookupObject<surfaceScalarField>(phiName_);

    tmp<surfaceScalarField> CoPtr(nullptr);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        CoPtr
            = mesh().surfaceInterpolation::deltaCoeffs()
            *(mag(phi)/mesh().magSf())
            *deltaT;
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            static_cast<const objectRegistry&>(mesh())
           .lookupObject<volScalarField>(rhoName_).oldTime();

        CoPtr =
            mesh().surfaceInterpolation::deltaCoeffs()
           *(mag(phi)/(fvc::interpolate(rho)*mesh().magSf()))
           *deltaT;
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);

        return tmp<surfaceScalarField>(nullptr);
    }

    if (mesh().time().timeIndex() < 2)
    {
        tmp<surfaceScalarField> rDt((CoPtr()/maxCo_)/deltaT);
        rDt->primitiveFieldRef() = max
        (
            scalar(1.0)/deltaT.value(),
            gMax(rDt->primitiveField())
        );

        return rDt;
    }
    else
    {
        return max(CoPtr()/maxCo_, scalar(1))/deltaT;
    }
}


template<class Type>
tmp<volScalarField> CoEulerDdtScheme<Type>::CorDeltaT
(
    const VolField<Type>& vf
) const
{
    const surfaceScalarField cofrDeltaT(CofrDeltaT(vf));

    tmp<volScalarField> tcorDeltaT
    (
        volScalarField::New
        (
            "CorDeltaT",
            cofrDeltaT.db(),
            mesh(),
            dimensionedScalar(cofrDeltaT.dimensions(), 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& corDeltaT = tcorDeltaT.ref();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, facei)
    {
        corDeltaT[owner[facei]] =
            max(corDeltaT[owner[facei]], cofrDeltaT[facei]);

        corDeltaT[neighbour[facei]] =
            max(corDeltaT[neighbour[facei]], cofrDeltaT[facei]);
    }

    const surfaceScalarField::Boundary& cofrDeltaTbf =
        cofrDeltaT.boundaryField();

    forAll(cofrDeltaTbf, patchi)
    {
        const fvsPatchScalarField& pcofrDeltaT = cofrDeltaTbf[patchi];
        const fvPatch& p = pcofrDeltaT.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pcofrDeltaT, patchFacei)
        {
            corDeltaT[faceCells[patchFacei]] = max
            (
                corDeltaT[faceCells[patchFacei]],
                pcofrDeltaT[patchFacei]
            );
        }
    }

    corDeltaT.correctBoundaryConditions();

    return tcorDeltaT;
}


template<class Type>
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::CofrDeltaT
(
    const VolField<Type>& vf
) const
{
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    const objectRegistry& obr(vf.db());
    const surfaceScalarField& phi =
        obr.lookupObject<surfaceScalarField>(phiName_);

    tmp<surfaceScalarField> CoPtr(nullptr);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        CoPtr
            = mesh().surfaceInterpolation::deltaCoeffs()
            *(mag(phi)/mesh().magSf())
            *deltaT;
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
           obr.lookupObject<volScalarField>(rhoName_).oldTime();

        CoPtr
            = mesh().surfaceInterpolation::deltaCoeffs()
            *(mag(phi)/(fvc::interpolate(rho)*mesh().magSf()))
            *deltaT;
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);

        return tmp<surfaceScalarField>(nullptr);
    }

    if (mesh().time().timeIndex() < 2)
    {
        tmp<surfaceScalarField> rDt((CoPtr()/maxCo_)/deltaT);
        rDt->primitiveFieldRef() = max
        (
            scalar(1.0)/deltaT.value(),
            gMax(rDt->primitiveField())
        );

        return rDt;
    }
    else
    {
        return max(CoPtr()/maxCo_, scalar(1))/deltaT;
    }
}


template<class Type>
tmp<VolField<Type>> CoEulerDdtScheme<Type>::fvcDdt(const dimensioned<Type>& dt)
{
    const word d2dt2name("ddt(" + dt.name() + ')');

    if (mesh().moving())
    {
        const volScalarField rDeltaT(CorDeltaT());

        tmp<VolField<Type>> tdtdt
        (
            VolField<Type>::New
            (
                d2dt2name,
                mesh(),
                dimensioned<Type>("0", dt.dimensions()/dimTime, Zero)
            )
        );

        tdtdt.ref().primitiveFieldRef() =
            rDeltaT.primitiveField()*dt.value()
           *(1.0 - mesh().Vsc0()/mesh().Vsc());

        return tdtdt;
    }
    else
    {
        return VolField<Type>::New
        (
            d2dt2name,
            mesh(),
            dimensioned<Type>("0", dt.dimensions()/dimTime, Zero),
            calculatedFvPatchField<Type>::typeName
        );
    }
}


template<class Type>
tmp<VolField<Type>> CoEulerDdtScheme<Type>::fvcDdt(const VolField<Type>& vf)
{
    const volScalarField rDeltaT(CorDeltaT(vf));

    const word ddtName("ddt(" + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT()*(vf() - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()),
            rDeltaT.boundaryField()
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
tmp<VolField<Type>> CoEulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT(vf));

    const word ddtName("ddt(" + rho.name() + ',' + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT()*rho*(vf() - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()),
            rDeltaT.boundaryField()*rho.value()
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
tmp<VolField<Type>> CoEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT(vf));
    const word ddtName("ddt(" + rho.name() + ',' + vf.name() + ')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT()
           *(
                rho()*vf()
              - rho.oldTime()()*vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.boundaryField()
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
tmp<VolField<Type>> CoEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT(vf));

    const word ddtName
    (
        "ddt(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')'
    );

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT()
           *(
                alpha()*rho()*vf()
              - alpha.oldTime()()*rho.oldTime()()
               *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.boundaryField()*
            (
                alpha.boundaryField()*rho.boundaryField()*vf.boundaryField()
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
               alpha*rho*vf - alpha.oldTime()*rho.oldTime()*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>> CoEulerDdtScheme<Type>::fvmDdt(const VolField<Type>& vf)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalarField rDeltaT(CorDeltaT(vf)().primitiveField());

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
tmp<fvMatrix<Type>>
CoEulerDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(CorDeltaT(vf)().primitiveField());

    fvm.diag() = rDeltaT*rho.value()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
CoEulerDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(CorDeltaT(vf)().primitiveField());

    fvm.diag() = rDeltaT*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
CoEulerDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(CorDeltaT(vf)().primitiveField());

    fvm.diag() =
        rDeltaT*alpha.primitiveField()*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf,
    const word interpolationName
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT(U)));

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
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
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const VolField<Type>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT(U)));

    fluxFieldType phiCorr
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
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf,
    const word interpolationName
)
{
    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT(U)));

        VolField<Type> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr
        (
            phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0, interpolationName)
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr)*rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == dimDensity*dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        return fvcDdtUfCorr(U, Uf, interpolationName);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);
        ::abort();
    }
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
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
            this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), phiCorr)*rDeltaT*phiCorr
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
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::meshPhi(const VolField<Type>&)
{
    tmp<surfaceScalarField> tmeshPhi
    (
        surfaceScalarField::New
        (
            "meshPhi",
            mesh(),
            dimensionedScalar(dimVolume/dimTime, 0)
        )
    );

    tmeshPhi->setOriented();
    return tmeshPhi;
}


template<class Type>
tmp<scalarField> CoEulerDdtScheme<Type>::meshPhi
(
    const VolField<Type>&,
    const label patchi
)
{
    return tmp<scalarField>
    (
        new scalarField(mesh().boundary()[patchi].size(), 0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

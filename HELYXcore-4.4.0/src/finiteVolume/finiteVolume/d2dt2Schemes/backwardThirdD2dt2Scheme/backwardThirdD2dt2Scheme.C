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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


#include "finiteVolume/d2dt2Schemes/backwardThirdD2dt2Scheme/backwardThirdD2dt2Scheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

template<class Type>
scalar backwardThirdD2dt2Scheme<Type>::deltaT_() const
{
    return mesh().time().deltaT().value();
}


template<class Type>
scalar backwardThirdD2dt2Scheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<class Type>
scalar backwardThirdD2dt2Scheme<Type>::deltaT0_
(
    const VolField<Type>& vf
) const
{
    if
    (
        vf.oldTime().timeIndex()
     == vf.oldTime().oldTime().timeIndex()
     //    vf.oldTime().oldTime().timeIndex()
     // == vf.oldTime().oldTime().oldTime().timeIndex()
    )
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
backwardThirdD2dt2Scheme<Type>::backwardThirdD2dt2Scheme(const fvMesh& mesh)
:
    d2dt2Scheme<Type>(mesh),
    ocCoeff_(new Function1Types::Constant<scalar>("ocCoeff", 1))
{}


template<class Type>
backwardThirdD2dt2Scheme<Type>::backwardThirdD2dt2Scheme(const fvMesh& mesh, Istream& is)
:
    d2dt2Scheme<Type>(mesh, is)
{
    token firstToken(is);

    if (firstToken.isNumber())
    {
        const scalar ocCoeff = firstToken.number();
        if (ocCoeff < 0 || ocCoeff > 1)
        {
            FatalIOErrorInFunction
            (
                is
            )   << "Blending coefficient = " << ocCoeff
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        ocCoeff_ = new Function1Types::Constant<scalar>
        (
            "ocCoeff",
            ocCoeff
        );
    }
    else
    {
        is.putBack(firstToken);
        dictionary dict(is);
        ocCoeff_ = Function1<scalar>::New("ocCoeff", dict);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> backwardThirdD2dt2Scheme<Type>::fvcD2dt2
(
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardThirdD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().topoChanged())
    {
        notImplemented(type() + ": not implemented for topo changes");
    }

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    // const scalar deltaT = deltaT_();
    // const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 11.0/6.0;
    const scalar coefft0 = 3.0;
    const scalar coefft00 = 3.0/2.0;
    const scalar coefft000 = 1.0/3.0;

    // return VolField<Type>::New
    // (
    //     "d2dt2(" + vf.name() + ')',
    //     rDeltaT*
    //     (
    //         coefft*backwardDdtScheme<Type>(mesh()).fvcDdt(vf)
    //       - coefft0*backwardDdtScheme<Type>(mesh()).fvcDdt(vf.oldTime())
    //       + coefft00*backwardDdtScheme<Type>
    //         (
    //             mesh()
    //         ).fvcDdt(vf.oldTime().oldTime())
    //     )
    // );

    return VolField<Type>::New
    (
        "d2dt2(" + vf.name() + ')',
        (1 - ocCoeff())*backwardD2dt2Scheme<Type>(mesh()).fvcD2dt2(vf)
      + ocCoeff()*sqr(rDeltaT)*
        (
            + sqr(coefft)*vf
            - 2.0*coefft*coefft0*vf.oldTime()
            + 2.0*coefft*coefft00*vf.oldTime().oldTime()
            - 2.0*coefft*coefft000*
                vf.oldTime().oldTime().oldTime()
            + sqr(coefft0)*vf.oldTime().oldTime()
            + 2.0*coefft0*coefft000*
                vf.oldTime().oldTime().oldTime().oldTime()
            - 2.0*coefft0*coefft00*
                vf.oldTime().oldTime().oldTime()
            + sqr(coefft00)*
                vf.oldTime().oldTime().oldTime().oldTime()
            - 2.0*coefft00*coefft000*
                vf.oldTime().oldTime().oldTime().oldTime().oldTime()
            + sqr(coefft000)*
                vf.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime()
        )
    );
}


template<class Type>
tmp<VolField<Type>> backwardThirdD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardThirdD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().topoChanged())
    {
        notImplemented(type() + ": not implemented for topo changes");
    }

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar coefft = 11.0/6.0;
    const scalar coefft0 = 3.0;
    const scalar coefft00 = 3.0/2.0;
    const scalar coefft000 = 1.0/3.0;

    return VolField<Type>::New
    (
        "d2dt2(" + vf.name() + ')',
        (1 - ocCoeff())*backwardD2dt2Scheme<Type>(mesh()).fvcD2dt2(vf)
      + ocCoeff()*sqr(rDeltaT)*
        (
            + sqr(coefft)*vf*rho
            - 2.0*coefft*coefft0*vf.oldTime()*rho.oldTime()
            + 2.0*coefft*coefft00*vf.oldTime().oldTime()*
                rho.oldTime().oldTime()
            - 2.0*coefft*coefft000*
                vf.oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime()
            + sqr(coefft0)*vf.oldTime().oldTime()*rho.oldTime().oldTime()
            + 2.0*coefft0*coefft000*
                vf.oldTime().oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime().oldTime()
            - 2.0*coefft0*coefft00*
                vf.oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime()
            + sqr(coefft00)*
                vf.oldTime().oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime()
            - 2.0*coefft00*coefft000*
                vf.oldTime().oldTime().oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime().oldTime()
            + sqr(coefft000)*
                vf.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime()*
                rho.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime()
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>> backwardThirdD2dt2Scheme<Type>::fvmD2dt2
(
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardThirdD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().topoChanged())
    {
        notImplemented(type() + ": not implemented for topo changes");
    }

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>(vf, vf.dimensions()*dimVol/dimTime/dimTime)
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/deltaT_();
    // const scalar deltaT = deltaT_();
    // const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 11.0/6.0;
    const scalar coefft0 = 3.0;
    const scalar coefft00 = 3.0/2.0;
    const scalar coefft000 = 1.0/3.0;

    fvm = (1 - ocCoeff())*backwardD2dt2Scheme<Type>(mesh()).fvmD2dt2(vf);

    fvm.diag() += ocCoeff()*(sqr(coefft)*rDeltaT
        *dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT))*mesh().V();

    fvm.source() += ocCoeff()*sqr(rDeltaT)*mesh().V()*
    (
      + 2.0*coefft*coefft0*vf.oldTime().primitiveField()
      - 2.0*coefft*coefft00*vf.oldTime().oldTime().primitiveField()
      + 2.0*coefft*coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft0)*vf.oldTime().oldTime().primitiveField()
      - 2.0*coefft0*coefft000*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
      + 2.0*coefft0*coefft00*vf.oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft00)*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
      + 2.0*coefft00*coefft000*vf.oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft000)*vf.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> backwardThirdD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardThirdD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().topoChanged())
    {
        notImplemented(type() + ": not implemented for topo changes");
    }

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/deltaT_();
    // const scalar deltaT = deltaT_();
    // const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 11.0/6.0;
    const scalar coefft0 = 3.0;
    const scalar coefft00 = 3.0/2.0;
    const scalar coefft000 = 1.0/3.0;

    fvm = (1 - ocCoeff())*backwardD2dt2Scheme<Type>(mesh()).fvmD2dt2(vf);

    fvm.diag() += ocCoeff()*(sqr(coefft)*rDeltaT
        *dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT))*rho*mesh().V();

    fvm.source() += ocCoeff()*sqr(rDeltaT)*rho*mesh().V()*
    (
      + 2.0*coefft*coefft0*vf.oldTime().primitiveField()
      - 2.0*coefft*coefft00*vf.oldTime().oldTime().primitiveField()
      + 2.0*coefft*coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft0)*vf.oldTime().oldTime().primitiveField()
      - 2.0*coefft0*coefft000*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
      + 2.0*coefft0*coefft00*vf.oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft00)*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
      + 2.0*coefft00*coefft000*vf.oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
      - sqr(coefft000)*vf.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> backwardThirdD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardThirdD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().topoChanged())
    {
        notImplemented(type() + ": not implemented for a topo changes");
    }

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/deltaT_();
    // const scalar deltaT = deltaT_();
    // const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 11.0/6.0;
    const scalar coefft0 = 3.0;
    const scalar coefft00 = 3.0/2.0;
    const scalar coefft000 = 1.0/3.0;

    fvm = (1 - ocCoeff())*backwardD2dt2Scheme<Type>(mesh()).fvmD2dt2(vf);

    fvm.diag() += ocCoeff()*(sqr(coefft)*rDeltaT
        *dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT))*rho*mesh().V();

    fvm.source() += ocCoeff()*sqr(rDeltaT)*mesh().V()*
    (
      + 2.0*coefft*coefft0*vf.oldTime().primitiveField()*rho.oldTime()
      - 2.0*coefft*coefft00*vf.oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime()
      + 2.0*coefft*coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime()
      - sqr(coefft0)*vf.oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime()
      - 2.0*coefft0*coefft000*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime().oldTime()
      + 2.0*coefft0*coefft00*vf.oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime()
      - sqr(coefft00)*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime().oldTime()
      + 2.0*coefft00*coefft000*vf.oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime().oldTime().oldTime()
      - sqr(coefft000)*vf.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime().primitiveField()
        *rho.oldTime().oldTime().oldTime().oldTime().oldTime().oldTime()
    );

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

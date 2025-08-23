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
    (c) Zeljko Tukovic
        Some additions by Philip Cardiff (solids4foam)
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/


#include "finiteVolume/d2dt2Schemes/backwardD2dt2Scheme/backwardD2dt2Scheme.H"
#include "finiteVolume/ddtSchemes/backwardDdtScheme/backwardDdtScheme.H"

#include "finiteVolume/d2dt2Schemes/EulerD2dt2Scheme/EulerD2dt2Scheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT_() const
{
    return mesh().time().deltaT().value();
}


template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT0_
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first time step
    // ZT
    // if (mesh().time().timeIndex() == 1)
    // {
    //     return EulerD2dt2Scheme<Type>(mesh()).fvcD2dt2(vf);
    // }

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    // ZT
    // const scalar coefft = 1.5;
    // const scalar coefft0 = 2.0;
    // const scalar coefft00 = 0.5;

    return VolField<Type>::New
    (
        "d2dt2(" + vf.name() + ')',
        rDeltaT*
        (
            coefft*backwardDdtScheme<Type>(mesh()).fvcDdt(vf)
          - coefft0*backwardDdtScheme<Type>(mesh()).fvcDdt(vf.oldTime())
          + coefft00*backwardDdtScheme<Type>
            (
                mesh()
            ).fvcDdt(vf.oldTime().oldTime())
        )
    );
}


template<class Type>
tmp<VolField<Type>> backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first time step
    // ZT
    // if (mesh().time().timeIndex() == 1)
    // {
    //     return EulerD2dt2Scheme<Type>(mesh()).fvcD2dt2(rho, vf);
    // }

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    // ZT
    // const scalar coefft = 1.5;
    // const scalar coefft0 = 2.0;
    // const scalar coefft00 = 0.5;

    return VolField<Type>::New
    (
        "d2dt2(" + vf.name() + ')',
        rDeltaT
       *(
            coefft*rho*backwardDdtScheme<Type>(mesh()).fvcDdt(vf)
          - coefft0*rho.oldTime()*backwardDdtScheme<Type>
            (
                mesh()
            ).fvcDdt(vf.oldTime())
          + coefft00*rho.oldTime().oldTime()*backwardDdtScheme<Type>
            (
                mesh()
            ).fvcDdt(vf.oldTime().oldTime())
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>> backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first time step
    // if (mesh().time().timeIndex() == 1)
    // {
    //     return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(vf);
    // }

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>(vf, vf.dimensions()*dimVol/dimTime/dimTime)
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/deltaT_();
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first time step
    // ZT
    // if (mesh().time().timeIndex() == 1)
    // {
    //     return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(rho, vf);
    // }

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
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*rho*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

    fvm.source() += rDeltaT*rho*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first time step
    // ZT
    // if (mesh().time().timeIndex() == 1)
    // {
    //     return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(rho, vf);
    // }

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
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*rho*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*rho.oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*rho.oldTime().oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

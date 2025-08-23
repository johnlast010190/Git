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
    (c) 2019-2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const VolField<Type>& vf,
    const word& name
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar(dimless, 1)
    );

    return fvm::laplacian(Gamma, vf, name);
}


template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const VolField<Type>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar(dimless, 1)
    );

    return fvm::laplacian
    (
        Gamma,
        vf,
        "laplacian(" + vf.name() + ')'
    );
}


template<class Type>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const VolField<Type>& vf,
    const word& name
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar(dimless, 1)
    );

    return fvm::bLaplacian(Gamma, vf, name);
}


template<class Type>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const VolField<Type>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar(dimless, 1)
    );

    return fvm::bLaplacian
    (
        Gamma,
        vf,
        "laplacian(" + vf.name() + ')'
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const zero&,
    const VolField<Type>& vf,
    const word& name
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const zero&,
    const VolField<Type>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const one&,
    const VolField<Type>& vf,
    const word& name
)
{
    return fvm::laplacian(vf, name);
}


template<class Type>
tmp<fvMatrix<Type>>
laplacian
(
    const one&,
    const VolField<Type>& vf
)
{
    return fvm::laplacian(vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    const SurfaceField<GType> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvm::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const VolField<Type>& vf
)
{
    const SurfaceField<GType> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvm::laplacian(Gamma, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Laplacian(fvm::laplacian(tgamma(), vf, name));
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvm::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> Laplacian(fvm::laplacian(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> tLaplacian = fvm::laplacian(tgamma(), vf, name);
    tgamma.clear();
    return tLaplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvm::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tGamma,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm(fvm::laplacian(tGamma(), vf));
    tGamma.clear();
    return tfvm;
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvm::bLaplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const tmp<SurfaceField<GType>>& tGamma,
    const VolField<Type>& vf
)
{
    tmp<fvBlockMatrix<Type>> tbfvm(fvm::bLaplacian(tGamma(), vf));
    tGamma.clear();
    return tbfvm;
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvm::bLaplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const tmp<VolField<GType>>& tGamma,
    const VolField<Type>& vf
)
{
    tmp<fvBlockMatrix<Type>> tbfvm(fvm::bLaplacian(tGamma(), vf));
    tGamma.clear();
    return tbfvm;
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().fvmBLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvBlockMatrix<Type>> tbLaplacian(fvm::bLaplacian(tgamma(), vf, name));
    tgamma.clear();
    return tbLaplacian;
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.db(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().fvmBLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp
<
    fvBlockMatrix<Type>
> bLaplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvBlockMatrix<Type>> tbLaplacian(fvm::bLaplacian(tgamma(), vf, name));
    tgamma.clear();
    return tbLaplacian;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2017-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/convectionSchemes/boundedConvectionScheme/boundedConvectionScheme.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const convectionScheme<Type>& boundedConvectionScheme<Type>::scheme() const
{
    return scheme_();
}


template<class Type>
tmp<SurfaceField<Type>>
boundedConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf
) const
{
    return scheme_().interpolate(phi, vf);
}


template<class Type>
tmp<SurfaceField<Type>>
boundedConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    return scheme_().flux(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvmDiv(faceFlux, vf)
      + fvm::SuSp(-fvc::surfaceIntegrate(absFlux), vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const volScalarField& vfMult,
    const VolField<Type>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvmDiv(faceFlux, vfMult, vf)
      + fvm::SuSp(-fvc::surfaceIntegrate(absFlux)*vfMult, vf);
}


template<class Type>
tmp<fvBlockMatrix<Type>>
boundedConvectionScheme<Type>::fvmBDiv
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    tmp<fvBlockMatrix<Type>> tbfvm
    (
        new fvBlockMatrix<Type>
        (
            const_cast<VolField<Type>&>(vf)
        )
    );

    tmp<fvMatrix<Type>> segConv = fvmDiv(faceFlux, vf);

    tbfvm->insertEquation(0, segConv.ref());

    return tbfvm;
}


template<class Type>
tmp<VolField<Type>>
boundedConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvcDiv(faceFlux, vf)
      - fvc::surfaceIntegrate(absFlux)*vf;
}


template<class Type>
tmp<surfaceScalarField> boundedConvectionScheme<Type>::absoluteFlux
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    if (faceFlux.mesh().moving() && faceFlux.dimensions() == dimVolume/dimTime)
    {
        return fvc::absolute(faceFlux, vf);
    }
    else
    {
        return tmp<surfaceScalarField>(faceFlux);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

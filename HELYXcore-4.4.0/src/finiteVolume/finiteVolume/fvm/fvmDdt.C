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
    (c) 2021 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const VolField<Type>& vf
)
{
    return fv::getDdtScheme(vf).ref().fvmDdt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvmDdt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    return fv::getDdtScheme(rho, vf).ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fv::getDdtScheme(rho, vf).ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const volScalarField& rho,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fv::getDdtScheme(alpha, rho, vf).ref().fvmDdt(alpha, rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const one&,
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const one&,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return ddt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
ddt
(
    const volScalarField& alpha,
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(alpha, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "finiteVolume/d2dt2Schemes/steadyStateD2dt2Scheme/steadyStateD2dt2Scheme.H"
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
tmp<VolField<Type>> steadyStateD2dt2Scheme<Type>::fvcD2dt2
(
    const VolField<Type>& vf
)
{
    return VolField<Type>::New
    (
        "d2dt2(" + vf.name() + ')',
        mesh(),
        dimensioned<Type>("0", vf.dimensions()/dimTime/dimTime, Zero)
    );
}


template<class Type>
tmp<VolField<Type>> steadyStateD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return VolField<Type>::New
    (
        "d2dt2(" + rho.name() + ',' + vf.name() + ')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            rho.dimensions()*vf.dimensions()/dimTime/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>> steadyStateD2dt2Scheme<Type>::fvmD2dt2
(
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>(vf, vf.dimensions()*dimVol/dimTime/dimTime)
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> steadyStateD2dt2Scheme<Type>::fvmD2dt2
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
            rho.dimensions()*vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>> steadyStateD2dt2Scheme<Type>::fvmD2dt2
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
            rho.dimensions()*vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

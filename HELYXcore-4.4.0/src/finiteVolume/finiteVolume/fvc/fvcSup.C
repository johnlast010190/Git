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

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> Su(const VolField<Type>& su, const VolField<Type>& vf)
{
    return su;
}


template<class Type>
tmp<VolField<Type>>
Su(const tmp<VolField<Type>>& tsu, const VolField<Type>& vf)
{
    return tsu;
}


template<class Type>
tmp<VolField<Type>> Sp(const volScalarField& sp, const VolField<Type>& vf)
{
    return sp*vf;
}


template<class Type>
tmp<VolField<Type>>
Sp(const tmp<volScalarField>& tsp, const VolField<Type>& vf)
{
    return tsp*vf;
}


template<class Type>
tmp<VolField<Type>>
Sp
(
    const dimensionedScalar& sp,
    const VolField<Type>& vf
)
{
    return sp*vf;
}


template<class Type>
tmp<VolField<Type>> SuSp(const volScalarField& sp, const VolField<Type>& vf)
{
    return sp*vf;
}


template<class Type>
tmp<VolField<Type>>
SuSp(const tmp<volScalarField>& tsp, const VolField<Type>& vf)
{
    return tsp*vf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

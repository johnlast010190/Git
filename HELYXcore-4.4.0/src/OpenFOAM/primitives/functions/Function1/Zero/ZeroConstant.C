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
    (c) 2017 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Zero/ZeroConstant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::ZeroConstant<Type>::ZeroConstant
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::ZeroConstant<Type>::~ZeroConstant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::ZeroConstant<Type>::value(const scalar x) const
{
    return pTraits<Type>::zero;
}


template<class Type>
Type Foam::Function1Types::ZeroConstant<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return pTraits<Type>::zero;
}


template<class Type>
Type Foam::Function1Types::ZeroConstant<Type>::derivative
(
    const scalar x
) const
{
    return pTraits<Type>::zero;
}


template<class Type>
void Foam::Function1Types::ZeroConstant<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::END_STATEMENT << nl;
}


// ************************************************************************* //

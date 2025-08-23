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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/NSRDS5/NSRDS5.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDS5<Type>::NSRDS5
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    Istream& is(dict.lookup(entryName));
    // Gobble type of function1
    word entryType(is);
    is >> coeffs_;
}


template<class Type>
Foam::Function1Types::NSRDS5<Type>::NSRDS5
(
    const word& entryName,
    const FixedList<Type, 4>& coeffs
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDS5<Type>::NSRDS5(const NSRDS5& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDS5<Type>::~NSRDS5()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::NSRDS5<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDS5<Type>::value(const scalar x) const
{
    // Evaluate the function
    // a_/pow(b_, 1 + pow(1 - T/c_, d_))
    // and return the result

    return
        cmptDivide
        (
            coeffs_[0],
            cmptPow
            (
                coeffs_[1],
                pTraits<Type>::one
                    + cmptPow
                    (
                        pTraits<Type>::one
                        - cmptDivide(pTraits<Type>::one*x, coeffs_[2]), coeffs_[3]
                    )
            )
        );
}


template<class Type>
Type Foam::Function1Types::NSRDS5<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDS5<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDS5<Type>::derivative(const scalar x) const
{
    // Evaluate the vectorised derivative of the function
    // f(x) = a_/pow(b_, 1 + pow(1 - x/c_, d_))
    // given by
    // f(x)*d_*log(b_)*((1 - x/c_)^(d_-1))/c_
    // and return the result

    return
        cmptMultiply
        (
            value(x),
            cmptDivide
            (
                cmptMultiply
                (
                    cmptLog
                    (
                        coeffs_[1]
                    ),
                    cmptMultiply
                    (
                        coeffs_[3],
                        cmptPow
                        (
                            pTraits<Type>::one
                          - cmptDivide
                            (
                                pTraits<Type>::one*x, coeffs_[2]
                            ),
                            coeffs_[3] - pTraits<Type>::one
                        )
                    )
                ), coeffs_[2]
            )
        );
}


template<class Type>
void Foam::Function1Types::NSRDS5<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    // Switch to ASCII for writing
    const bool isBinary(os.format() == IOstream::BINARY);
    isBinary ? os.format(IOstream::ASCII) : os.format(IOstream::ASCII);
    os  << nl << indent << coeffs_
        << token::END_STATEMENT << nl;
    isBinary ? os.format(IOstream::BINARY) : os.format(IOstream::ASCII);
}


// ************************************************************************* //

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
    (c) 2020-2024 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/NonUniformTable/NonUniformTable1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NonUniformTable<Type>::NonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    Function1<Type>(name),
    low_(GREAT),
    high_(-GREAT),
    values_(),
    delta_(GREAT)
{
    if (values_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Table " << nl
            << "    " << name << nl
            << "    has less than 2 entries."
            << exit(FatalIOError);
    }
    else
    {
        low_ = values_.first().first();
        high_ = values_.last().first();

        for (label i = 1; i<values_.size(); i++)
        {
            delta_ = min(delta_, values_[i].first() - values_[i - 1].first());
        }

        delta_ *= 0.9;

        jumpTable_.setSize((high_ - low_)/delta_ + 1);

        label i = 0;
        forAll(jumpTable_, j)
        {
            const scalar x = low_ + j*delta_;

            if (x > values_[i + 1].first())
            {
                i++;
            }

            jumpTable_[j] = i;
        }
    }
}


template<class Type>
Foam::Function1Types::NonUniformTable<Type>::NonUniformTable
(
    const NonUniformTable<Type>& nut
)
:
    Function1<Type>(nut),
    low_(nut.low_),
    high_(nut.high_),
    values_(nut.values_),
    delta_(nut.delta_),
    jumpTable_(nut.jumpTable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::NonUniformTable<Type>::value
(
    scalar x
) const
{
    const label i = index(x);
    const scalar xi = values_[i].first();
    const scalar lambda = (x - xi)/(values_[i + 1].first() - xi);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


template<class Type>
Foam::scalar Foam::Function1Types::NonUniformTable<Type>::integrate
(
    const scalar x
) const
{
    const label i = index(x);
    const scalar xi = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar dx = x - xi;
    const scalar lambda = dx/(values()[i + 1].first() - xi);

    return (fi + 0.5*lambda*(values()[i + 1].second() - fi))*dx;
}


template<class Type>
Foam::scalar Foam::Function1Types::NonUniformTable<Type>::integrateByX
(
    const scalar x
) const
{
    const label i = index(x);
    const scalar xi = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar gradf =
        (values()[i + 1].second() - fi)/(values()[i + 1].first() - xi);

    return ((fi - gradf*xi)*log(x/xi) + gradf*(x - xi));
}


template<class Type>
Foam::scalar Foam::Function1Types::NonUniformTable<Type>::integrateByX
(
    const scalar x1,
    const scalar x2
) const
{
    return integrateByX(x2) - integrateByX(x1);
}


template<class Type>
Type Foam::Function1Types::NonUniformTable<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NonUniformTable<Type>::derivative
(
    scalar x
) const
{
    const label i = index(x);

    return
        (values_[i + 1].second() - values_[i].second())
       /(values_[i + 1].first() - values_[i].first());
}


template<class Type>
void Foam::Function1Types::NonUniformTable<Type>::writeData(Ostream& os) const
{
    os.writeEntry("values", values_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::NonUniformTable<Type>::operator=
(
    const NonUniformTable<Type>& nut
)
{
    low_ = nut.low_;
    high_ = nut.high_;
    values_ = nut.values_;
    delta_ = nut.delta_;
    jumpTable_ = nut.jumpTable_;
}


// ************************************************************************* //

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
    (c) 2019-2010 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/uniformTable1/uniformTable1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::UniformTable<Type>::UniformTable
(
    const word& name,
    const dictionary& dict
)
:
    Function1<Type>(name),
    dictName_(dict.name()),
    low_(dict.lookup<scalar>("low")),
    high_(dict.lookup<scalar>("high")),
    values_(dict.lookup("values"))
{
    if (values_.size() < 2)
    {
        FatalErrorInFunction
            << "Table " << nl
            << "    " << dictName_ << nl
            << "    has less than 2 entries."
            << exit(FatalError);
    }
    else
    {
        delta_ = (high_ - low_)/(values_.size() - 1);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::UniformTable<Type>::value(scalar x) const
{
    const scalar nd = (x - low_)/delta_;
    const label i = nd;

    if (nd < 0 || i > values_.size() - 2)
    {
        FatalErrorInFunction
            << x << " out of range "
            << low_ << " to " << high_ << nl
            << "    of table " << dictName_
            << exit(FatalError);
    }

    const scalar xi = low_ + i*delta_;
    const scalar lambda = (x - xi)/delta_;

    return values_[i] + lambda*(values_[i + 1] - values_[i]);
}


template<class Type>
Type Foam::Function1Types::UniformTable<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
void Foam::Function1Types::UniformTable<Type>::writeData(Ostream& os) const
{
    os.writeEntry("low", low_);
    os.writeEntry("high", high_);
    os.writeEntry("values", values_);
}


// ************************************************************************* //

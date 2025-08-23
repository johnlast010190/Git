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
    (c) 2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::label Foam::Function1Types::NonUniformTable<Type>::index
(
    scalar x
) const
{
    if (x < low_ || x > high_)
    {
        FatalErrorInFunction
            << x << " out of range "
            << low_ << " to " << high_ << nl
            << "    of NonUniformTable " << this->name_
            << exit(FatalError);
    }

    const scalar nd = (x - low_)/delta_;
    const label j = nd;

    label i = jumpTable_[j];

    if (i < values_.size() - 1 && x > values_[i + 1].first())
    {
        i++;
    }

    return i;
}


// ************************************************************************* //

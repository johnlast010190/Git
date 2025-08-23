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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pointProcConstraintData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pointProcConstraintData<Type>::pointProcConstraintData
(
    const bool& pp,
    const bool& pbp,
    const bool& fixesValue,
    const Type& fixedValue
)
:
    procPoint_(pp),
    procBoundaryPoint_(pbp),
    fixesValue_(fixesValue),
    fixedValue_(fixedValue)
{
}


template<class Type>
Foam::pointProcConstraintData<Type>::pointProcConstraintData()
:
    procPoint_(false),
    procBoundaryPoint_(false),
    fixesValue_(false),
    fixedValue_(Zero)
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
//
/*
template<class Type>
void Foam::pointProcConstraintData<Type>::operator=
(
    const Foam::pointProcConstraintData<Type>& rhs
) const
{
    procPoint() == rhs.procPoint();
    procBoundaryPoint() == rhs.procBoundaryPoint();
    fixesValue() = rhs.fixesValue();
    fixedValue() = rhs.fixedValue();
}
*/

template<class Type>
bool Foam::pointProcConstraintData<Type>::operator==
(
    const Foam::pointProcConstraintData<Type>& rhs
) const
{
    return
    (
        procPoint() == rhs.procPoint()
     && procBoundaryPoint() == rhs.procBoundaryPoint()
     && fixesValue() == rhs.fixesValue()
     && fixedValue() == rhs.fixedValue()
    );
}


template<class Type>
bool Foam::pointProcConstraintData<Type>::operator!=
(
    const Foam::pointProcConstraintData<Type>& rhs
) const
{
    return !(*this == rhs);
}


template<class Type>
bool Foam::pointProcConstraintData<Type>::operator>
(
    const Foam::pointProcConstraintData<Type>& rhs
) const
{
    if (!fixesValue() && rhs.fixesValue())
    {
        return false;
    }
    else
    {
        if (fixesValue() && rhs.fixesValue())
        {
            //- Check if both fix value and select the one with the
            //  larger fixedValue. It should happen rarely
            if (mag(fixedValue()) > mag(rhs.fixedValue()))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return true;
        }
    }
}


template<class Type>
bool Foam::pointProcConstraintData<Type>::operator<
(
    const Foam::pointProcConstraintData<Type>& rhs
) const
{
    return !(*this < rhs);
}


template<class Type>
Foam::pointProcConstraintData<Type> Foam::max
(
    const Foam::pointProcConstraintData<Type>& a1,
    const Foam::pointProcConstraintData<Type>& a2
)
{
    return (a1 > a2) ? a1 : a2;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * *


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, pointProcConstraintData<Type>& iD)
{
    is  >> iD.procPoint()
        >> iD.procBoundaryPoint()
        >> iD.fixesValue()
        >> iD.fixedValue();

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>(Foam::Istream&, "
        "Foam::pointProcConstraintData&)"
    );

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pointProcConstraintData<Type>& iD
)
{
    os  << iD.procPoint() << token::SPACE
        << iD.procBoundaryPoint() << token::SPACE
        << iD.fixesValue() << token::SPACE
        << iD.fixedValue() << token::SPACE
        << endl;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::pointProcConstraintData&)"
    );

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

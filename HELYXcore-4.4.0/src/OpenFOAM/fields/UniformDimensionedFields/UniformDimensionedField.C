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
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/UniformDimensionedFields/UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensioned<Type>& dt
)
:
    regIOobject(io),
    dimensioned<Type>(dt)
{
    // Read value
    readHeaderOk(IOstream::BINARY, typeName);
}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dictionary& dict
)
:
    regIOobject(io),
    dimensioned<Type>(regIOobject::name(), dimless, pTraits<Type>::zero)
{
    this->dimensions().reset(dict.lookup("dimensions"));
    this->value() = dict.lookup("value");
}



template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const UniformDimensionedField<Type>& rdt
)
:
    regIOobject(rdt),
    dimensioned<Type>(rdt)
{}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io
)
:
    regIOobject(io),
    dimensioned<Type>(regIOobject::name(), dimless, Zero)
{
    // For if MUST_READ_IF_MODIFIED
    addWatch();

    // Read unless NO_READ
    readHeaderOk(IOstream::BINARY, typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::~UniformDimensionedField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::UniformDimensionedField<Type>::readData(Istream& is)
{
    dictionary dict(is);
    scalar multiplier;
    this->dimensions().read(dict.lookup("dimensions"), multiplier);
    this->value() = dict.lookup<Type>("value");
    this->value() *= multiplier;

    return is.good();
}


template<class Type>
bool Foam::UniformDimensionedField<Type>::writeData(Ostream& os) const
{
    scalar multiplier;
    os.writeKeyword("dimensions");
    this->dimensions().write(os, multiplier) << token::END_STATEMENT << nl;
    os.writeEntry("value", (this->value()/multiplier));
    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const UniformDimensionedField<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const dimensioned<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


// ************************************************************************* //

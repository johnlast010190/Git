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
    (c) 2017-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "thermophysicalProperties/thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalProperties, 0);
    defineRunTimeSelectionTable(thermophysicalProperties,);
    defineRunTimeSelectionTable(thermophysicalProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalProperties::thermophysicalProperties(scalar W)
:
    W_(W)
{}


Foam::thermophysicalProperties::thermophysicalProperties(const dictionary& dict)
:
    W_
    (
        dict.found("W")
      ? dict.lookup<scalar>("W")
      : dict.lookup<scalar>("molWeight")
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    const auto ctor =
        ctorTableLookup
        (
            "thermophysicalProperties type",
            ConstructorTable_(),
            name
        );
    return autoPtr<thermophysicalProperties>(ctor());
}


Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    const word& thermophysicalPropertiesTypeName = dict.dictName();

    const auto ctor =
        ctorTableLookup
        (
            "thermophysicalProperties type",
            dictionaryConstructorTable_(),
            thermophysicalPropertiesTypeName
        );
    return autoPtr<thermophysicalProperties>(ctor(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermophysicalProperties::readIfPresent(const dictionary& dict)
{
    dict.found("W")
  ? dict.readIfPresent("W", W_)
  : dict.readIfPresent("molWeight", W_);
}


void Foam::thermophysicalProperties::write(Ostream& os) const
{
    os.writeEntry("W", W_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const thermophysicalProperties& l)
{
    l.write(os);
    return os;
}


// ************************************************************************* //

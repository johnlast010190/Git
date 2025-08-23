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
    (c) 2015-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOobjects/IOdictionary/localIOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localIOdictionary::localIOdictionary(const IOobject& io)
:
    IOdictionary(io, typeName)
{
    readHeaderOk(IOstream::ASCII, typeName);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    IOdictionary(io, typeName)
{
    if (!readHeaderOk(IOstream::ASCII, typeName))
    {
        dictionary::operator=(dict);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    const IOobject& io,
    const word& actualType
)
:
    IOdictionary(io, actualType)
{
    readHeaderOk(IOstream::ASCII, actualType);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


// ************************************************************************* //

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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Andrade/Andrade.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(Andrade);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::Andrade::Andrade
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<scalar>(entryName)
{
    Istream& is(dict.lookup(entryName));
    // Gobble type of function1
    word entryType(is);
    is >> coeffs_;
}


Foam::Function1Types::Andrade::Andrade(const Andrade& se)
:
    Function1<scalar>(se),
    coeffs_(se.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::Andrade::~Andrade()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::Andrade::value(const scalar T) const
{
    return exp
    (
        coeffs_[0]
      + coeffs_[1]*T
      + coeffs_[2]*sqr(T)
      + coeffs_[3]/(coeffs_[4] + T)
    );
}


void Foam::Function1Types::Andrade::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);

    // Switch to ASCII for writing
    const bool isBinary(os.format() == IOstream::BINARY);
    os.format(IOstream::ASCII);

    os  << nl << indent << coeffs_ << token::END_STATEMENT << nl;

    isBinary ? os.format(IOstream::BINARY) : os.format(IOstream::ASCII);
}


// ************************************************************************* //

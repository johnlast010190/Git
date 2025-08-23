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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "polyTopoChange/refinementData/refinementData.H"

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::refinementData& wDist
)
{
    if (os.format() == IOstream::ASCII)
    {
        os << wDist.refinementCount_ << token::SPACE << wDist.count_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&wDist.refinementCount_),
            sizeof(refinementData)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Foam::Istream& is, Foam::refinementData& wDist)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> wDist.refinementCount_ >> wDist.count_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&wDist.refinementCount_),
            sizeof(refinementData)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //

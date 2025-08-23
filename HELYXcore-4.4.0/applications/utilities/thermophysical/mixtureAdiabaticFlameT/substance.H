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

Class
    Foam::substance

Description

SourceFiles
    substanceI.H
    substance.C
    substanceIO.C

\*---------------------------------------------------------------------------*/

#ifndef substance_H
#define substance_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/word/word.H"
#include "db/IOstreams/IOstreams/Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class substance Declaration
\*---------------------------------------------------------------------------*/

class substance
{
    // Private data

        word name_;
        scalar volFrac_;


public:

    // Constructors

        //- Construct null
        substance()
        {}


    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }

            scalar volFrac() const
            {
                return volFrac_;
            }

            bool operator==(const substance& s) const
            {
                return name_ == s.name_ && volFrac_ == s.volFrac_;
            }

            bool operator!=(const substance& s) const
            {
                return !operator==(s);
            }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, substance& s)
        {
            is  >> s.name_ >> s.volFrac_;
            return is;
        }

        friend Ostream& operator<<(Ostream& os, const substance& s)
        {
            os  << s.name_ << token::SPACE << s.volFrac_;
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

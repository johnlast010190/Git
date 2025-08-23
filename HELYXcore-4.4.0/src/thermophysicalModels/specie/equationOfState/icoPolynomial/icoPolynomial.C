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
    (c) 2015 Engys Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/icoPolynomial/icoPolynomial.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize>::icoPolynomial
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    Specie(obr, dict),
    rhoCoeffs_
    (
        dict.subDict("equationOfState").lookup
        (
            "rhoCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    rhoMin_
    (
        dict.subDict("equationOfState").lookupOrDefault<scalar>("rhoMin", 0.0)
    ),
    rhoMax_
    (
        dict.subDict("equationOfState").lookupOrDefault<scalar>("rhoMax", GREAT)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
void Foam::icoPolynomial<Specie, PolySize>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add
    (
        word("rhoCoeffs<" + Foam::name(PolySize) + '>'),
        rhoCoeffs_
    );

    if (rhoMin_ > 0)
    {
        dict.add(word("rhoMin"), rhoMin_);
    }
    if (rhoMax_ < GREAT)
    {
        dict.add(word("rhoMax"), rhoMax_);
    }

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const icoPolynomial<Specie, PolySize>& ip
)
{
    ip.write(os);
    return os;
}


// ************************************************************************* //

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

#include "correlationFunction/correlationFunction.H"
#include "db/IOstreams/IOstreams.H"

template<class Type>
bool Foam::correlationFunction<Type>::writeAveraged(Ostream& os) const
{
    Field<scalar> averageCF(averaged());

    forAll(averageCF, v)
    {
        os  << v*sampleInterval()
            << token::SPACE
            << averageCF[v]
            << nl;
    }

    return os.good();
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const correlationFunction<Type>& cF
)
{
    os  << cF.duration()
        << nl << cF.sampleInterval()
        << nl << cF.averagingInterval()
        << nl << cF.sampleSteps()
        << nl << cF.tZeroBuffers()
        << nl << static_cast<const bufferedAccumulator<scalar>&>(cF);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

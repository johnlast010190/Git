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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "algorithms/PointEdgeWave/PointData.H"

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class DataType>
Foam::Ostream& Foam::operator<<(Ostream& os, const PointData<DataType>& pd)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << static_cast<const pointEdgePoint&>(pd)
            << token::SPACE << pd.data();
    }
    else
    {
        return os
            << static_cast<const pointEdgePoint&>(pd)
            << pd.data();
    }
}


template<class DataType>
Foam::Istream& Foam::operator>>(Istream& is, PointData<DataType>& pd)
{
    return is >> static_cast<pointEdgePoint&>(pd) >> pd.data_;
}


// ************************************************************************* //

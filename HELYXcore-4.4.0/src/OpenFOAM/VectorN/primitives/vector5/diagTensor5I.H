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
    (c) 2017-2024 Engys Ltd.

Type
    diagTensor5

Description
    diagonalTensorN of 5 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //


//- Inner-product between a diagonal tensor and a vector
inline vector
operator&(const diagTensor5& st, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and diagTensor5 "
        << "vector has 3 components while diagTensor5 has 5"
        << exit(FatalError);

    return result;
}

//- Inner-product between a diagonal tensor and a vector
inline vector
operator&(const vector& v, const diagTensor5& st)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and diagTensor5 "
        << "vector has 3 components while diagTensor5 has 5"
        << exit(FatalError);


    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

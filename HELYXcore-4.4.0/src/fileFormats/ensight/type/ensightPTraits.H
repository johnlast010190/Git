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
    (c) 2013-2015 OpenFOAM Foundation

Class
    Foam::ensightPTraits

Description
    Conversion of OpenFOAM pTraits into the Ensight equivalent

\*---------------------------------------------------------------------------*/

#ifndef ensightPTraits_H
#define ensightPTraits_H

#include "primitives/pTraits/pTraits.H"
#include "fields/Fields/fieldTypes.H"
#include "primitives/direction/direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightPTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class ensightPTraits
{
public:

    // Static data members

        static const char* const typeName;

        //- Ordering table: return OpenFOAM component given Ensight component
        //  This is used for the symmTensor ordering: Ensight requires
        //  xx yy zz xy xz yz
        static const direction componentOrder[];
};

template<>
const char* const ensightPTraits<scalar>::typeName;

template<>
const direction ensightPTraits<scalar>::componentOrder[];

template<>
const char* const ensightPTraits<vector>::typeName;

template<>
const direction ensightPTraits<vector>::componentOrder[];

template<>
const char* const ensightPTraits<sphericalTensor>::typeName;

template<>
const direction ensightPTraits<sphericalTensor>::componentOrder[];

template<>
const char* const ensightPTraits<symmTensor>::typeName;

template<>
const direction ensightPTraits<symmTensor>::componentOrder[];

template<>
const char* const ensightPTraits<tensor>::typeName;

template<>
const direction ensightPTraits<tensor>::componentOrder[];


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

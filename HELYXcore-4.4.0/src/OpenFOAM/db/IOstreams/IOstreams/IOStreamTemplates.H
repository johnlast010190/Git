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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef HELYXCORE_IOSTREAMTEMPLATES_HPP
#define HELYXCORE_IOSTREAMTEMPLATES_HPP

#include <tuple>

// A place for IO operators that should be defined universally, and should be always inlineable
// for performance reasons. Generally pertains to STL entities.


namespace Foam
{


/*---------------------------------------------------------------------------*\
          Tuple IO
\*---------------------------------------------------------------------------*/


template<size_t... Is, typename... Ts>
void writeTuple(
    std::index_sequence<Is...>,
    Ostream& os,
    const std::tuple<Ts...>& tuple
) {
    ((os << std::get<Is>(tuple)), ...);
}

template<size_t... Is, typename... Ts>
void readTuple(
    std::index_sequence<Is...>,
    Istream& os,
    std::tuple<Ts...>& tuple
) {
    ((os >> std::get<Is>(tuple)), ...);
}


template<typename... Ts>
Ostream& operator<<(Ostream& os, const std::tuple<Ts...>& tuple) {
    writeTuple(std::index_sequence_for<Ts...>{}, os, tuple);
    return os;
}

template<typename... Ts>
Istream& operator>>(Istream& is, std::tuple<Ts...>& tuple) {
    readTuple(std::index_sequence_for<Ts...>{}, is, tuple);
    return is;
}

}


#endif //HELYXCORE_IOSTREAMTEMPLATES_HPP

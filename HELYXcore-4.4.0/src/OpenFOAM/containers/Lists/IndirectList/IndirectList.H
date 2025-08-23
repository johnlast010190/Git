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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::IndirectList

Description
    A List with indirect addressing.

See also
    Foam::UIndirectList for a version without any allocation for the
    addressing.

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef IndirectList_H
#define IndirectList_H

#include "containers/Lists/UIndirectList/UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IndirectListAddressing Declaration
\*---------------------------------------------------------------------------*/

//- A helper class for storing addresses.
class IndirectListAddressing
{
    // Private data

        //- Storage for the list addressing
        List<label> addressing_;


protected:

    // Constructors

        //- Construct by copying the addressing array
        explicit inline IndirectListAddressing(const labelUList& addr);

        //- Move constructor by transferring addressing array
        explicit inline IndirectListAddressing(List<label>&& addr);

        //- Disallow default bitwise copy construction
        IndirectListAddressing(const IndirectListAddressing&) = delete;


    // Member Functions

        // Access

            //- Return the list addressing
            inline const List<label>& addressing() const;

        // Edit

            //- Reset addressing
            inline void resetAddressing(const labelUList&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IndirectListAddressing&) = delete;
};


/*---------------------------------------------------------------------------*\
                        Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
:
    private IndirectListAddressing,
    public  UIndirectList<T>
{
public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline IndirectList(const UList<T>&, const labelUList&);

        //- Move constructor given the complete list and moves the addressing
        inline IndirectList(const UList<T>&, List<label>&&);

        //- Copy constructor
        inline IndirectList(const IndirectList<T>&);

        //- Construct from UIndirectList
        explicit inline IndirectList(const UIndirectList<T>&);


    // Member Functions

        // Access

            //- Return the list addressing
            using UIndirectList<T>::addressing;


        // Edit

            //- Reset addressing
            using IndirectListAddressing::resetAddressing;

            //- Copy the pointer held by the given IndirectList.
            inline void shallowCopy(const IndirectList<T>&);


        // Member Operators

            //- Disallow default assignment operator
            void operator=(const IndirectList<T>&) = delete;

            //- Disallow assignment from UIndirectList
            void operator=(const UIndirectList<T>&) = delete;

            //- Assignment operator
            using UIndirectList<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/Lists/IndirectList/IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2012 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

Class
    Foam::IOMRFSourceList

Description
    List of MRF sources with IO functionality.

SourceFiles
    IOMRFSourceList.C

\*---------------------------------------------------------------------------*/

#ifndef IOMRFSourceList_H
#define IOMRFSourceList_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "MRFSourceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace  fv
{

/*---------------------------------------------------------------------------*\
                        Class IOMRFSourceList Declaration
\*---------------------------------------------------------------------------*/

class IOMRFSourceList
:
    public IOdictionary,
    public MRFSourceList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const objectRegistry& obr) const;


public:

    // Constructors

        //- Construct from registry
        IOMRFSourceList(const objectRegistry& obr);

        //- Disallow default bitwise copy construction
        IOMRFSourceList(const IOMRFSourceList&) = delete;


        //- Destructor
        virtual ~IOMRFSourceList()
        {}


    // Member Functions

        //- Read dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IOMRFSourceList&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

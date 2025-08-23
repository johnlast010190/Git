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
    (c) 2010-2015 Engys Ltd.

Class
    Foam::offsetLduInterface

Description
    An abstract base class to be shared by offsetGAMGInterface and
    offsetLduInterface

SourceFiles
    offsetLduInterfaceBase.C

\*---------------------------------------------------------------------------*/

#ifndef offsetLduInterfaceBase_H
#define offsetLduInterfaceBase_H

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "fields/Fields/primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class offsetLduInterfaceBase Declaration
\*---------------------------------------------------------------------------*/

class offsetLduInterfaceBase
{

public:

    //- Runtime type information
    TypeName("offsetLduInterfaceBase");

    // Constructors

    //- Destructor
    virtual ~offsetLduInterfaceBase();


    // Member Functions

        // Access

            //- Return controlled interface
            virtual const lduInterface& slaveInterface() const = 0;

//            virtual const label internalOffset() const = 0;

//            virtual const label internalSize() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2024 Engys Ltd.

Class
    Foam::lduInterface

Description
    An abstract base class for implicitly-coupled interfaces
    e.g. processor and cyclic patches.

SourceFiles
    lduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterface_H
#define lduInterface_H

#include "fields/Fields/labelField/labelField.H"
#include "db/typeInfo/typeInfo.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                     Class lduInterface Declaration
\*---------------------------------------------------------------------------*/

class lduInterface
{
public:

    //- Runtime type information
    TypeName("lduInterface");


    // Constructors

        //- Construct null
        lduInterface()
        {}

        //- Disallow default bitwise copy construction
        lduInterface(const lduInterface&) = delete;


    //- Destructor
    virtual ~lduInterface();


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const = 0;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {}

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const = 0;

            //- Transfer and return internal field adjacent to the interface
            //  Can be used when a derived class needs to supply interfaces
            //  that should behave differently when being run in segregated
            //  vs region-coupled mode
            virtual tmp<labelField> regionCoupledInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {
                return internalFieldTransfer(commsType, iF);
            }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lduInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

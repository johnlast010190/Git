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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::offsetLduInterface

Description

SourceFiles
    offsetLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef offsetLduInterface_H
#define offsetLduInterface_H

#include "fields/Fields/labelField/labelField.H"
#include "db/typeInfo/typeInfo.H"
#include "db/IOstreams/Pstreams/Pstream.H"

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "offsetLduInterfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                     Class offsetLduInterface Declaration
\*---------------------------------------------------------------------------*/

class offsetLduInterface
:
    public lduInterface,
    virtual public offsetLduInterfaceBase
{
    // Private Member Variables

        //- Base interface that we are controlling
        const lduInterface& slaveInterface_;

        //- Offset into internal cell array
        const label internalOffset_;

        //- Size of internal cell array
        const label internalSize_;

        //- Offset into neighbouring cell array
        const label foreignOffset_;

        //- Size of neighbouring cell array
        const label foreignSize_;

        //- Persistent storage for return of ref
        mutable labelList faceCells_;


public:

    //- Runtime type information
    TypeName("offset");


    // Constructors

        //- Construct taking patches from the field
        offsetLduInterface
        (
            const lduInterface& interface,
            const label internalOffset,
            const label internalSize,
            const label foreignOffset,
            const label foreignSize
        );

        //- Disallow default bitwise copy construction
        offsetLduInterface(const offsetLduInterface&) = delete;


    //- Destructor
    virtual ~offsetLduInterface();


    // Member Functions

        // Access

            //- Return underlying interface we are controlling
            virtual const lduInterface& slaveInterface() const
            {
                return slaveInterface_;
            }

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const offsetLduInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

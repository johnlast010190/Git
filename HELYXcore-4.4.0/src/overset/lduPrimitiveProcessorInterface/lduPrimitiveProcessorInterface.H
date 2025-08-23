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
    (c) 2016 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::lduPrimitiveProcessorInterface

Description
    Concrete implementation of processor interface. Used to temporarily store
    settings.

SourceFiles
    lduPrimitiveProcessorInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveProcessorInterface_H
#define lduPrimitiveProcessorInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class lduPrimitiveProcessorInterface Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveProcessorInterface
:
    public lduInterface,
    public processorLduInterface
{
    // Private data

        const labelList faceCells_;

        //- My processor rank in communicator
        label myProcNo_;

        //- Neighbouring processor rank in communicator
        label neighbProcNo_;

        //- Transformation tensor
        tensorField forwardT_;

        //- Message tag used for sending
        int tag_;

        //- Communicator to use for parallel communication
        const label comm_;


public:

    //- Runtime type information
    TypeName("primitiveProcessor");


    // Constructors

        //- Construct from components
        lduPrimitiveProcessorInterface
        (
            const labelUList& faceCells,
            const label myProcNo,
            const label neighbProcNo,
            const tensorField& forwardT,
            const int tag,
            const label comm = UPstream::worldComm
        );

        //- Disallow default bitwise copy construction
        lduPrimitiveProcessorInterface
        (
            const lduPrimitiveProcessorInterface&
        ) = delete;


    //- Destructor
    virtual ~lduPrimitiveProcessorInterface()
    {}


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return faceCells_;
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise neighbour field transfer
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


        //- Processor interface functions

            //- Return communicator used for sending
            virtual label comm() const
            {
                return comm_;
            }

            //- Return processor number (rank in communicator)
            virtual int myProcNo() const
            {
                return myProcNo_;
            }

            //- Return neigbour processor number (rank in communicator)
            virtual int neighbProcNo() const
            {
                return neighbProcNo_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return dummy null transform
            virtual const transformer& transform() const
            {
                return transformer::null;
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return tag_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveProcessorInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

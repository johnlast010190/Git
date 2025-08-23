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

Class
    Foam::nonConformalDiscreteMixingGAMGInterface

Description
    GAMG agglomerated interface for non-conformal discrete mixing patches.

SourceFiles
    nonConformalDiscreteMixingGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingGAMGInterface_H
#define nonConformalDiscreteMixingGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class discreteMixingIntersection;

/*---------------------------------------------------------------------------*\
                   Class nonConformalLduInterface Declaration
\*---------------------------------------------------------------------------*/

class nonConformalLduInterface
{
public:

    //- Runtime type information
    ClassNameNoDebug("nonConformalLduInterface");

    //- Destructor
    virtual ~nonConformalLduInterface(){}

    //- Return neighbour patch ID
    virtual label nbrPatchID() const = 0;

    //- Does this side own the patch ?
    virtual bool owner() const = 0;

    //- Return neighbour patch interface
    virtual const nonConformalLduInterface& nbrPatch() const = 0;

    //- Access the list of intersection engines
    virtual const discreteMixingIntersection& intersection() const = 0;
};


/*---------------------------------------------------------------------------*\
           Class nonConformalDiscreteMixingGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class nonConformalDiscreteMixingGAMGInterface
:
    public GAMGInterface,
    virtual public nonConformalLduInterface
{
    // Private Data

        //- Reference to the nonConformalLduInterface from which this is
        //  agglomerated
        const nonConformalLduInterface& fineInterface_;

        //- Non-conformal discrete mixing interface
        autoPtr<discreteMixingIntersection> intersectionsPtr_;

        //- Offset from the patch ID in the original matrix to that in the
        //  super-matrix
        const label patchIDOffset_;


public:

    //- Runtime type information
    TypeName("nonConformalDiscreteMixing");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        nonConformalDiscreteMixingGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Disallow default bitwise copy construction
        nonConformalDiscreteMixingGAMGInterface
        (
            const nonConformalDiscreteMixingGAMGInterface&
        ) = delete;


    //- Destructor
    virtual ~nonConformalDiscreteMixingGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const override;


        // Non-conformal interface functions

            //- Return neighbour patch ID
            virtual label nbrPatchID() const override
            {
                return (fineInterface_.nbrPatchID() + patchIDOffset_);
            }

            //- Does this side own the patch ?
            virtual bool owner() const override
            {
                return fineInterface_.owner();
            }

            //- Return neighbour patch
            virtual const
            nonConformalDiscreteMixingGAMGInterface& nbrPatch() const override
            {
                return refCast<const nonConformalDiscreteMixingGAMGInterface>
                (
                    coarseInterfaces_[nbrPatchID()]
                );
            }

            //- Return the intersection engine
            virtual const
            discreteMixingIntersection& intersection() const override
            {
                return intersectionsPtr_();
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const override
            {
                NotImplemented;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonConformalDiscreteMixingGAMGInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

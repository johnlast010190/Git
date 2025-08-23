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
    (c) 2017 OpenCFD Ltd.

Class
    Foam::oversetGAMGInterface

Description
    GAMG agglomerated overset interface.

SourceFiles
    oversetGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetGAMGInterface_H
#define oversetGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "oversetPolyPatch/oversetLduInterface.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oversetGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetGAMGInterface
:
    public GAMGInterface,
    virtual public oversetLduInterface
{
    // Private data

        //- Reference for the oversetLduInterface from which this is
        //  agglomerated
        const oversetLduInterface& fineOversetInterface_;

        //- Temporary copy of the fine level restrict map. Cleared upon
        //  calculating stencils below
        mutable labelList restrictMap_;

        labelListList stencil_;

        autoPtr<distributionMap> mapPtr_;

        List<scalarList> cellInterpolationWeights_;

        scalarField normalisation_;

        labelList interpolationCells_;

        scalarList cellInterpolationWeight_;


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        oversetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Disallow default bitwise copy construction
        oversetGAMGInterface(const oversetGAMGInterface&) = delete;


    //- Destructor
    virtual ~oversetGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Overset interface functions

            //- Name of lduInterface
            virtual const word& name() const
            {
                return fineOversetInterface_.name();
            }

            //- Am I the master interface
            virtual bool master() const
            {
                return fineOversetInterface_.master();
            }

            //- GAMG restriction (fine-to-coarse)
            virtual const labelList& restrictMap() const
            {
                return restrictMap_;
            }

            //- Donor stencil
            virtual const labelListList& stencil() const
            {
                return stencil_;
            }

            //- Map for obtaining data in stencil order
            virtual const distributionMap& cellInterpolationMap() const
            {
                return mapPtr_();
            }

            //- Weights in stencil order
            virtual const List<scalarList>& cellInterpolationWeights() const
            {
                return cellInterpolationWeights_;
            }

            //- Normalisation ot matrix; for explicit contributions
            virtual const scalarField& normalisation() const
            {
                return normalisation_;
            }

            //- Acceptor cells
            virtual const labelList& interpolationCells() const
            {
                return interpolationCells_;
            }

            //- Underrelaxation for acceptor cells
            virtual const scalarList& cellInterpolationWeight() const
            {
                return cellInterpolationWeight_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                //TBD.
                NotImplemented;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const oversetGAMGInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

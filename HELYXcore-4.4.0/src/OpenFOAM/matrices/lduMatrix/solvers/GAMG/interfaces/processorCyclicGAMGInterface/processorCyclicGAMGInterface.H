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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::processorCyclicGAMGInterface

Description
    GAMG agglomerated processor interface.

SourceFiles
    processorCyclicGAMGInterface.C
    processorCyclicGAMGInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicGAMGInterface_H
#define processorCyclicGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/processorGAMGInterface/processorGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorCyclicGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicGAMGInterface
:
    public processorGAMGInterface
{
public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        processorCyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from Istream
        processorCyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        processorCyclicGAMGInterface
        (
            const processorCyclicGAMGInterface&
        ) = delete;


    // Destructor

        virtual ~processorCyclicGAMGInterface();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const processorCyclicGAMGInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

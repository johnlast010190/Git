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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::processorCyclicGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorCyclicGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicGAMGInterfaceField_H
#define processorCyclicGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/processorGAMGInterfaceField/processorGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class processorCyclicGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicGAMGInterfaceField
:
    public processorGAMGInterfaceField
{
public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        processorCyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );

        //- Construct from GAMG interface and fine level interface field
        processorCyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        processorCyclicGAMGInterfaceField
        (
            const processorCyclicGAMGInterfaceField&
        ) = delete;


    // Destructor

        virtual ~processorCyclicGAMGInterfaceField();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const processorCyclicGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

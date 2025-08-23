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
    Foam::GAMGInterfaceField

Description
     Abstract base class for GAMG agglomerated interface fields.

SourceFiles
    GAMGInterfaceField.C
    GAMGInterfaceFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGInterfaceField_H
#define GAMGInterfaceField_H

#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.H"
#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class GAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class GAMGInterfaceField
:
    public lduInterfaceField
{
    // Private data

        //- Local reference cast into the interface
        const GAMGInterface& interface_;


public:

    //- Runtime type information
    TypeName("GAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterfaceField,
            lduInterfaceField,
            (
                const GAMGInterface& GAMGCp,
                const lduInterfaceField& fineInterface
            ),
            (GAMGCp, fineInterface)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterfaceField,
            lduInterface,
            (
                const GAMGInterface& GAMGCp,
                const int rank
            ),
            (GAMGCp, rank)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterfaceField> New
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterfaceField> New
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        GAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField&
        )
        :
            lduInterfaceField(GAMGCp),
            interface_(GAMGCp)
        {}

        //- Construct from GAMG interface and fine level interface field
        GAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        )
        :
            lduInterfaceField(GAMGCp),
            interface_(GAMGCp)
        {}

        //- Disallow default bitwise copy construction
        GAMGInterfaceField(const GAMGInterfaceField&) = delete;


    // Member Functions

        //- Return interface
        const GAMGInterface& interface() const
        {
            return interface_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

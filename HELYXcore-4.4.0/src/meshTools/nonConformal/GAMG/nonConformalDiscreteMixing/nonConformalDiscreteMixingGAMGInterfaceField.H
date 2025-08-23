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
    Foam::nonConformalDiscreteMixingGAMGInterfaceField

Description
    GAMG agglomerated interface field for non-conformal discrete mixing patches.

SourceFiles
    nonConformalDiscreteMixingGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingGAMGInterfaceField_H
#define nonConformalDiscreteMixingGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class nonConformalDiscreteMixingGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class nonConformalDiscreteMixingGAMGInterfaceField
:
    public GAMGInterfaceField
{
    // Private Data

        //- Local reference cast into the non-conformal interface
        const nonConformalDiscreteMixingGAMGInterface& nonConformalInterface_;


public:

    //- Runtime type information
    TypeName("nonConformalDiscreteMixing");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        nonConformalDiscreteMixingGAMGInterfaceField
        (
            const GAMGInterface& GAMGNonConformalP,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and component rank
        nonConformalDiscreteMixingGAMGInterfaceField
        (
            const GAMGInterface& GAMGNonConformalP,
            const int rank
        );

        //- Disallow default bitwise copy construction
        nonConformalDiscreteMixingGAMGInterfaceField
        (
            const nonConformalDiscreteMixingGAMGInterfaceField&
        ) = delete;


    //- Destructor
    virtual ~nonConformalDiscreteMixingGAMGInterfaceField();


    // Member Functions

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const bool add,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonConformalDiscreteMixingGAMGInterfaceField&) =
            delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

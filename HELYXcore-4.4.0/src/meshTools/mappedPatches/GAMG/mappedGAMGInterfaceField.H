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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::mappedGAMGInterfaceField

Description
    GAMG agglomerated region coupled interface field.

SourceFiles
    mappedGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedGAMGInterfaceField_H
#define mappedGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "mappedPatches/GAMG/mappedGAMGInterface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mappedGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class mappedGAMGInterfaceField
:
    public GAMGInterfaceField
{
    // Private data

        //- Local reference cast into the region coupled interface
        const mappedGAMGInterface& regionCoupledGAMGInterface_;

        //- Agglomerated explicit correction (for radiation)
        scalarField explicitNbrCorr_;


public:

    //- Runtime type information
    TypeName("mappedPatch");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        mappedGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        mappedGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        mappedGAMGInterfaceField(const mappedGAMGInterfaceField&) = delete;


    //- Destructor
    virtual ~mappedGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return regionCoupledGAMGInterface_.size();
            }

            virtual tmp<scalarField> explicitNbrCorr() const
            {
                return tmp<scalarField>(new scalarField(explicitNbrCorr_));
            }


    //- Interface matrix update
            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const bool add,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mappedGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

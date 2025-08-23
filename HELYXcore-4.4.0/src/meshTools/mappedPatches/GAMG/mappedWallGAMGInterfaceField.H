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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::mappedWallGAMGInterfaceField

Description
    GAMG agglomerated region coupled wall interface field.

SourceFiles
    mappedWallGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedWallGAMGInterfaceField_H
#define mappedWallGAMGInterfaceField_H

#include "mappedPatches/GAMG/mappedGAMGInterfaceField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mappedGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class mappedWallGAMGInterfaceField
:
    virtual public mappedGAMGInterfaceField
{
public:

    //- Runtime type information
    TypeName("mappedWall");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        mappedWallGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        mappedWallGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        mappedWallGAMGInterfaceField
        (
            const mappedWallGAMGInterfaceField&
        ) = delete;


    //- Destructor
    virtual ~mappedWallGAMGInterfaceField();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mappedWallGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

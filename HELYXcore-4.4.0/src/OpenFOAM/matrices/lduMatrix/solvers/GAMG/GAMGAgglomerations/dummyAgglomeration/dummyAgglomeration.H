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
    (c) 2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::dummyAgglomeration

Description
    Agglomerate without combining cells. Used for testing.

SourceFiles
    dummyAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef dummyAgglomeration_H
#define dummyAgglomeration_H

#include "matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/GAMGAgglomeration/GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dummyAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class dummyAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        //- Preset number of levels
        label nLevels_;


public:

    //- Runtime type information
    TypeName("dummy");


    // Constructors

        //- Construct given mesh and controls
        dummyAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Disallow default bitwise copy construction
        dummyAgglomeration(const dummyAgglomeration&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dummyAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

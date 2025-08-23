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
    (c) 2016 Engys Ltd.

Class
    Foam::tHelyxSegregatedEulerEulerLESState

Description
    interFoam state identifiers and specialised assembly functions

SourceFiles
    tHelyxSegregatedEulerEulerLESState.C

\*---------------------------------------------------------------------------*/

#ifndef tHelyxSegregatedEulerEulerLESState_H
#define tHelyxSegregatedEulerEulerLESState_H

#include "stateFunctions/tHelyxSegregatedEulerEulerRAS/tHelyxSegregatedEulerEulerRASState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tHelyxSegregatedEulerEulerLESState Declaration
\*---------------------------------------------------------------------------*/

class tHelyxSegregatedEulerEulerLESState
:
    public tHelyxSegregatedEulerEulerRASState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tHelyxSegregatedEulerEulerLESState(const tHelyxSegregatedEulerEulerLESState&);

        //- Disallow default bitwise assignment
        void operator=(const tHelyxSegregatedEulerEulerLESState&);


public:

    //- Runtime type information
        TypeName("tHelyxSegregatedEulerEulerLES");

    // Constructors

        //- Construct from components
        tHelyxSegregatedEulerEulerLESState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tHelyxSegregatedEulerEulerLESState()
    {}


    // Member Functions

        // Access

            // type indicators

            // turbulence (RANS/LES/none)
            inline virtual turbulenceType turbulence() const
            {
                return tuLES;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

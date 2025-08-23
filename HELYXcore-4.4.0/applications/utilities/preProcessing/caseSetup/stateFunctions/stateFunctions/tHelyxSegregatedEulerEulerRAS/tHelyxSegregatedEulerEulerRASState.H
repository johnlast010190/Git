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
    Foam::tHelyxSegregatedEulerEulerRASState

Description
    interFoam state identifiers and specialised assembly functions

SourceFiles
    tHelyxSegregatedEulerEulerRASState.C

\*---------------------------------------------------------------------------*/

#ifndef tHelyxSegregatedEulerEulerRASState_H
#define tHelyxSegregatedEulerEulerRASState_H

#include "stateComponents/EulerEuler/EulerEulerState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
               Class tHelyxSegregatedEulerEulerRASState Declaration
\*---------------------------------------------------------------------------*/

class tHelyxSegregatedEulerEulerRASState
:
    public EulerEulerState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tHelyxSegregatedEulerEulerRASState(const tHelyxSegregatedEulerEulerRASState&);

        //- Disallow default bitwise assignment
        void operator=(const tHelyxSegregatedEulerEulerRASState&);


public:

    //- Runtime type information
        TypeName("tHelyxSegregatedEulerEulerRAS");

    // Constructors

        //- Construct from components
        tHelyxSegregatedEulerEulerRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tHelyxSegregatedEulerEulerRASState()
    {}


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline virtual regionType region() const
            {
                return rtFluid;
            }

            // time (transient/steady)
            inline virtual timeType time() const
            {
                return ttTrans;
            }

            // turbulence (RANS/LES/none)
            inline virtual turbulenceType turbulence() const
            {
                return tuRAS;
            }

            // compressibility (compressible/incompressible/none)
            inline virtual compressibilityType compressibility() const
            {
                return ctComp;
            }

            //- USF vs legacy solver state
            virtual bool usf() const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

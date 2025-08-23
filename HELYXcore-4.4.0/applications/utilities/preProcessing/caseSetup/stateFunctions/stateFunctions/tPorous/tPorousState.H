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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::tPorousState

Description

SourceFiles
    tPorousState.C

\*---------------------------------------------------------------------------*/

#ifndef tPorousState_H
#define tPorousState_H

#include "stateComponents/porousSolid/porousSolidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tPorousState Declaration
\*---------------------------------------------------------------------------*/

class tPorousState
:
    public porousSolidState
{
public:

    //- Runtime type information
    TypeName("tPorous");

    // Constructors

        //- Construct from components
        tPorousState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        tPorousState(const tPorousState&) = delete;


    //- Destructor
    virtual ~tPorousState()
    {}


    // Member Functions

        //- Region (fluid/solid/surface/none)
        inline virtual regionType region() const override
        {
            return rtSolid;
        }

        //- Time (transient/steady)
        inline virtual timeType time() const override
        {
            return ttSteady;
        }

        //- Turbulence (RANS/LES/none)
        inline virtual turbulenceType turbulence() const override
        {
            return tuNone;
        }

        //- Compressibility (compressible/incompressible/none)
        inline virtual compressibilityType compressibility() const override
        {
            return ctNone;
        }

        //- USF vs legacy solver state
        virtual bool usf() const override
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tPorousState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2016-2024 Engys Ltd.

Class
    Foam::chtTransientSolidState

Description
    Transient solid state identifiers and specialised assembly functions

SourceFiles
    chtTransientSolidState.C

\*---------------------------------------------------------------------------*/

#ifndef chtTransientSolidState_H
#define chtTransientSolidState_H

#include "stateComponents/solid/solidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class chtTransientSolidState Declaration
\*---------------------------------------------------------------------------*/

class chtTransientSolidState
:
    public solidState
{
public:

    //- Runtime type information
    TypeName("chtTransientSolid");

    // Constructors

        //- Construct from components
        chtTransientSolidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        chtTransientSolidState(const chtTransientSolidState&) = delete;


    //- Destructor
    virtual ~chtTransientSolidState()
    {}


    // Member Functions

        //- Time (transient/steady)
        inline virtual timeType time() const
        {
            return ttTrans;
        }

        //- Turbulence (RANS/LES/none)
        inline virtual turbulenceType turbulence() const
        {
            return tuNone;
        }

        //- Compressibility (compressible/incompressible/none)
        inline virtual compressibilityType compressibility() const
        {
            return ctNone;
        }

        //- USF vs legacy solver state
        virtual bool usf() const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const chtTransientSolidState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

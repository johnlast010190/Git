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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::tHelyxCoupledRASState

Description
    helyxCoupled state identifiers and specialised assembly functions

SourceFiles
    tHelyxCoupledRASState.C

\*---------------------------------------------------------------------------*/

#ifndef tHelyxCoupledRASState_H
#define tHelyxCoupledRASState_H

#include "stateFunctions/helyxCoupled/helyxCoupledState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tHelyxCoupledRASState Declaration
\*---------------------------------------------------------------------------*/

class tHelyxCoupledRASState
:
    public helyxCoupledState
{
public:

    //- Runtime type information
    TypeName("tHelyxCoupledRAS");

    // Constructors

        //- Construct from components
        tHelyxCoupledRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        tHelyxCoupledRASState(const tHelyxCoupledRASState&) = delete;


    //- Destructor
    virtual ~tHelyxCoupledRASState()
    {}


    // Member Functions

        //- Time (transient/steady)
        inline virtual timeType time() const
        {
            return ttTrans;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tHelyxCoupledRASState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::buoyantMulticomponentRASState

Description

SourceFiles
    buoyantMulticomponentRASState.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantMulticomponentRASState_H
#define buoyantMulticomponentRASState_H

#include "stateComponents/multicomponentFluid/multicomponentFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class buoyantMulticomponentRASState Declaration
\*---------------------------------------------------------------------------*/

class buoyantMulticomponentRASState
:
    public multicomponentFluidState
{
public:

    //- Runtime type information
    TypeName("buoyantMulticomponentRAS");

    // Constructors

        //- Construct from components
        buoyantMulticomponentRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        buoyantMulticomponentRASState
        (
            const buoyantMulticomponentRASState&
        ) = delete;


    //- Destructor
    virtual ~buoyantMulticomponentRASState()
    {}


    // Member Functions

        //- Region (fluid/solid/surface/none)
        inline virtual regionType region() const
        {
            return rtFluid;
        }

        //- Time (transient/steady)
        inline virtual timeType time() const
        {
            return ttTrans;
        }

        //- Turbulence (RANS/LES/none)
        inline virtual turbulenceType turbulence() const
        {
            return tuRAS;
        }

        //- Compressibility (compressible/incompressible/none)
        inline virtual compressibilityType compressibility() const
        {
            return ctComp;
        }

        //- Final manipulations
        virtual void finalise();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantMulticomponentRASState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

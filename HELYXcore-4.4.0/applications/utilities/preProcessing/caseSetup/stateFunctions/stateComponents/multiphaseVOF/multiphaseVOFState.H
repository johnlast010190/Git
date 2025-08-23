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
    Foam::multiphaseVOFState

Description
    Abstract, multiPhaseVOF state transport/thermo assembly

SourceFiles
    multiphaseVOFState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseVOFState_H
#define multiphaseVOFState_H

#include "stateComponents/multiphaseFluid/multiphaseFluidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseVOFState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseVOFState
:
    public multiphaseFluidState
{
    // Private Functions

        //- initialise two phase (multiphaseVOF)
        void initialiseTwoPhase();

        //- initialise multiphase (multiphasemultiphaseVOF)
        void initialiseMultiPhase();


public:

    //- Runtime type information
    TypeName("multiphaseVOF");

    // Constructors

        //- Construct from components
        multiphaseVOFState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        multiphaseVOFState(const multiphaseVOFState&) = delete;


    //- Destructor
    virtual ~multiphaseVOFState()
    {}


    // Member Functions

        //- Assemble field definitions
        virtual void initialise();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseVOFState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

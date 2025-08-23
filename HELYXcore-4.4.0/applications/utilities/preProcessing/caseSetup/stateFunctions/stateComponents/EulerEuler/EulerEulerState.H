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
    Foam::EulerEulerState

Description
    Abstract, Euler-Euler state for USF - will eventually replace
    multiphaseEuler

SourceFiles
    EulerEulerState.C

\*---------------------------------------------------------------------------*/

#ifndef EulerEulerState_H
#define EulerEulerState_H

#include "stateComponents/multiphaseMulticomponent/multiphaseMulticomponentState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class EulerEulerState Declaration
\*---------------------------------------------------------------------------*/

class EulerEulerState
:
    public multiphaseMulticomponentState
{

protected:

    // Protected Data

    // Protected member functions


private:

    // Private Data


    // Private Functions

        //- Disallow default bitwise copy construct
        EulerEulerState(const EulerEulerState&);

        //- Disallow default bitwise assignment
        void operator=(const EulerEulerState&);



public:

    //- Runtime type information
        TypeName("EulerEuler");

    // Constructors

        //- Construct from components
        EulerEulerState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~EulerEulerState()
    {}


    // Member Functions

        // Access

            //- Provide the names of the phases
            virtual wordList phaseNames() const override
            {
                return materialNames();
            }


        // Check

        // Edit

            //- assemble field definitions
            virtual void initialise() override;

            //- assemble field definitions
            virtual void finalise() override;

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

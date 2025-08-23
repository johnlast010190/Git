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
    Foam::VOFState

Description
    Abstract, VOF state for USF - will eventually replace multiphaseVOF

SourceFiles
    VOFState.C

\*---------------------------------------------------------------------------*/

#ifndef VOFState_H
#define VOFState_H

#include "stateComponents/multiphaseMulticomponent/multiphaseMulticomponentState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class VOFState Declaration
\*---------------------------------------------------------------------------*/

class VOFState
:
    public multiphaseMulticomponentState
{
public:

    //- Runtime type information
    TypeName("VOF");

    // Constructors

        //- Construct from components
        VOFState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        VOFState(const VOFState&) = delete;


    //- Destructor
    virtual ~VOFState()
    {}


    // Member Functions

        //- Provide the names of the phases
        virtual wordList phaseNames() const override
        {
            return materialNames();
        }

        //- Assemble field definitions
        virtual void initialise() override;

        //- Assemble field definitions
        virtual void finalise() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VOFState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::multiphaseEuler

Description
    Abstract, multiphaseEuler state transport/thermo assembly

SourceFiles
    multiphaseEuler.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseEuler_H
#define multiphaseEuler_H

#include "stateComponents/multiphaseFluid/multiphaseFluidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseEuler Declaration
\*---------------------------------------------------------------------------*/

class multiphaseEulerState
:
    public multiphaseFluidState
{
    // Private Functions

        //- Scrape diameter model for inclusion into phase properties
        dictionary eulerPhaseProperties
        (
            const word& materialName,
            const dictionary&
        ) const;

        //- Initialise multiphase, might extend to twophaseEuler in future
        void initialiseMultiPhase();


public:

    //- Runtime type information
    TypeName("multiphaseEuler");

    // Constructors

        //- Construct from components
        multiphaseEulerState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        multiphaseEulerState(const multiphaseEulerState&) = delete;


    //- Destructor
    virtual ~multiphaseEulerState()
    {}


    // Member Functions

        //- Assemble field definitions
        virtual void initialise();

        //- Remove default U entry from fieldMaps
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseEulerState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

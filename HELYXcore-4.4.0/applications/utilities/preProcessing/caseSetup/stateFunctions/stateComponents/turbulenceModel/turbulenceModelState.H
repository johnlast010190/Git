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
    Foam::turbulenceModelState

Description
    turbulenceModel state identifiers and specialised assembly functions

SourceFiles
    turbulenceModelState.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceModelState_H
#define turbulenceModelState_H

#include "stateComponents/region/regionState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class turbulenceModelState Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModelState
:
    public regionState
{
protected:

    // Protected Data

        // Names of turbulence fields for each phase
        List<wordList> turbulenceFields_;

        //- Whether turbulence modelling is per-phase (as in Euler-Euler)
        const bool phasicTurbulence_;

    // Protected member functions

        virtual primitiveEntry divScheme(const word&, const dictionary&) const;

        static void addTurbulenceSpecificSchemes(const dictionary& turbDict, dictionary& fvSchemes);


private:

    // Private Functions

        //- Add solver and relaxation entries
        void fvSolutionInject
        (
            const dictionary& fieldDict,
            const dictionary& modelDict,
            dictionary& fvSolDict,
            const word& tfname
        );


public:

    //- Runtime type information
    TypeName("turbulenceModel");

    // Constructors

        //- Construct from components
        turbulenceModelState
        (
            const word& region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            const word& meshName = word::null,
            const bool phasicTurbulence = false
        );

        //- Disallow default bitwise copy construction
        turbulenceModelState(const turbulenceModelState&) = delete;


    //- Destructor
    ~turbulenceModelState() override = default;


    // Member Functions

        // region (fluid/solid/surface/none)
        inline regionType region() const override
        {
            return rtFluid;
        }

        //- assemble field definitions
        void initialise() override;

        //- Correct stateDict
        void correct() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModelState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

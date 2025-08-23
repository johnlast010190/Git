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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::noBlending

Description

SourceFiles
    noBlending.C

\*---------------------------------------------------------------------------*/

#ifndef noBlending_H
#define noBlending_H

#include "../eulerianBlendingMethod/eulerianBlendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianBlendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class noBlending Declaration
\*---------------------------------------------------------------------------*/

class noBlending
:
    public eulerianBlendingMethod
{
    // Private data

        //- Name of the continuous phase
        const word continuousPhase_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and two phases
        noBlending
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~noBlending();


    // Member Functions

        //- Factor for primary phase
        virtual tmp<volScalarField> f1
        (
            const eulerianPhaseModel& phase1,
            const eulerianPhaseModel& phase2
        ) const;

        //- Factor for secondary phase
        virtual tmp<volScalarField> f2
        (
            const eulerianPhaseModel& phase1,
            const eulerianPhaseModel& phase2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianBlendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::orderedEulerianPhasePair

Description

SourceFiles
    orderedEulerianPhasePair.C

\*---------------------------------------------------------------------------*/

#ifndef orderedPhasePair_H
#define orderedPhasePair_H

#include "../eulerianPhasePair/eulerianPhasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class orderedEulerianPhasePair Declaration
\*---------------------------------------------------------------------------*/

class orderedEulerianPhasePair
:
    public eulerianPhasePair
{

public:

    // Constructors

        //- Construct from two phases and gravity
        orderedEulerianPhasePair
        (
            const eulerianPhaseModel& dispersed,
            const eulerianPhaseModel& continuous
        );


    //- Destructor
    virtual ~orderedEulerianPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const eulerianPhaseModel& dispersed() const;

        //- Continuous phase
        virtual const eulerianPhaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

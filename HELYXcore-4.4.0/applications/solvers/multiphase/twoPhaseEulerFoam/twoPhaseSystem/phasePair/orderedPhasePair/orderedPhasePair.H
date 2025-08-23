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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::orderedPhasePair

Description

SourceFiles
    orderedPhasePair.C

\*---------------------------------------------------------------------------*/

#ifndef orderedPhasePair_H
#define orderedPhasePair_H

#include "phasePair/phasePair/phasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

/*---------------------------------------------------------------------------*\
                         Class orderedPhasePair Declaration
\*---------------------------------------------------------------------------*/

class orderedPhasePair
:
    public phasePair
{
    // Private data

        //- Aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;


public:

    // Constructors

        //- Construct from two phases, gravity, surface tension and aspect
        //  ratio tables
        orderedPhasePair
        (
            const phaseModel& dispersed,
            const phaseModel& continuous,
            const dimensionedVector& g,
            const scalarTable& sigmaTable,
            const dictTable& aspectRatioTable
        );


    //- Destructor
    virtual ~orderedPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const phaseModel& dispersed() const;

        //- Continuous phase
        virtual const phaseModel& continuous() const;

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

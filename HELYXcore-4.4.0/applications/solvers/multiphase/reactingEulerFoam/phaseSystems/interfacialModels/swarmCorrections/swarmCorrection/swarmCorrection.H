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
    Foam::swarmCorrection

Description

SourceFiles
    swarmCorrection.C
    newSwarmCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef swarmCorrection_H
#define swarmCorrection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class swarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class swarmCorrection
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("swarmCorrection");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        swarmCorrection,
        dictionary,
        (
            const dictionary& dict,
            const phasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and a phase pair
        swarmCorrection
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~swarmCorrection();


    // Selectors

        static autoPtr<swarmCorrection> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

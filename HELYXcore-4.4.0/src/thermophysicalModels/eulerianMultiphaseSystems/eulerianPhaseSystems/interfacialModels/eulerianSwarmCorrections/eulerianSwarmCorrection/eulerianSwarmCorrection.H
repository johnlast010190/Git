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
    Foam::eulerianSwarmCorrection

Description

SourceFiles
    eulerianSwarmCorrection.C
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

class eulerianPhasePair;

/*---------------------------------------------------------------------------*\
                           Class eulerianSwarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class eulerianSwarmCorrection
{
protected:

    // Protected data

        //- Phase pair
        const eulerianPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("swarmCorrectionModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        eulerianSwarmCorrection,
        dictionary,
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and a phase pair
        eulerianSwarmCorrection
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~eulerianSwarmCorrection();


    // Selectors

        static autoPtr<eulerianSwarmCorrection> New
        (
            const word& modelType,
            const dictionary& dict,
            const eulerianPhasePair& pair
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

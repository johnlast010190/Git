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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::departureFrequencyModel

Description
    Base class for bubble departure frequency models

SourceFiles
    departureFrequencyModel.C
    newdepartureFrequencyModel.C

\*---------------------------------------------------------------------------*/

#ifndef departureFrequencyModel_H
#define departureFrequencyModel_H

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

#include "phaseModel/phaseModel/phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class departureFrequencyModel Declaration
\*---------------------------------------------------------------------------*/

class departureFrequencyModel
{

public:

    //- Runtime type information
    TypeName("departureFrequencyModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        departureFrequencyModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        departureFrequencyModel();

        //- Disallow default bitwise copy construction
        departureFrequencyModel(const departureFrequencyModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<departureFrequencyModel> New(const dictionary& dict);


    //- Destructor
    virtual ~departureFrequencyModel();


    // Member Functions

        //- Calculate and return the bubble departure frequency
        virtual tmp<scalarField> fDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& dDep
        ) const = 0;

        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const departureFrequencyModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

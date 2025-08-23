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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::nucleationSiteModel

Description
    Base class for nucleation site density models

SourceFiles
    nucleationSiteModel.C
    newnucleationSiteModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationSiteModel_H
#define nucleationSiteModel_H

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
                         Class nucleationSiteModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationSiteModel
{

public:

    //- Runtime type information
    TypeName("nucleationSiteModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationSiteModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        nucleationSiteModel();

        //- Disallow default bitwise copy construction
        nucleationSiteModel(const nucleationSiteModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<nucleationSiteModel> New(const dictionary& dict);


    //- Destructor
    virtual ~nucleationSiteModel();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> N
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const = 0;

        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nucleationSiteModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

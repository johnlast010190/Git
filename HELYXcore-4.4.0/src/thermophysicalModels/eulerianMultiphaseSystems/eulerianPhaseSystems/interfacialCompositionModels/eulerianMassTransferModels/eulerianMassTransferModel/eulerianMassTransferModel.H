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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::eulerianMassTransferModel

Description

SourceFiles
    eulerianMassTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef massTransferModel_H
#define massTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

namespace Foam
{

class eulerianPhasePair;

/*---------------------------------------------------------------------------*\
                        Class eulerianMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class eulerianMassTransferModel
{
protected:

    // Protected data

        //- Phase pair
        const eulerianPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("massTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianMassTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const eulerianPhasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        //- Construct from a dictionary and a phase pair
        eulerianMassTransferModel
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~eulerianMassTransferModel();


    // Selectors

        static autoPtr<eulerianMassTransferModel> New
        (
            const word& modelType,
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    // Member Functions

        //- The implicit mass transfer coefficient
        //  Note: this has had the species mass diffusivity factored out
        virtual tmp<volScalarField> K() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

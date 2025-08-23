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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::kineticTheoryModels::granularPressureModel

SourceFiles
    granularPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef granularPressureModel_H
#define granularPressureModel_H

#include "db/dictionary/dictionary.H"
#include "fields/volFields/volFields.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class granularPressureModel Declaration
\*---------------------------------------------------------------------------*/

class granularPressureModel
{
protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("granularPressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        granularPressureModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        granularPressureModel(const dictionary& dict);

        //- Disallow default bitwise copy construction
        granularPressureModel(const granularPressureModel&) = delete;


    // Selectors

        static autoPtr<granularPressureModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~granularPressureModel();


    // Member Functions

        //- Granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeff
        (
            const volScalarField& alpha1,
            const volScalarField& g0,
            const volScalarField& rho1,
            const dimensionedScalar& e
        ) const = 0;

        //- Derivative of the granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeffPrime
        (
            const volScalarField& alpha1,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const volScalarField& rho1,
            const dimensionedScalar& e
        ) const = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const granularPressureModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

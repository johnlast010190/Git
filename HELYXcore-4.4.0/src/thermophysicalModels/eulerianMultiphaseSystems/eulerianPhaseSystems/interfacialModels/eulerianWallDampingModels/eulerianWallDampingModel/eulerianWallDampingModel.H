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
    Foam::eulerianWallDampingModel

Description

SourceFiles
    eulerianWallDampingModel.C
    newWallDampingModel.C

\*---------------------------------------------------------------------------*/

#ifndef wallDampingModel_H
#define wallDampingModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "../../eulerianWallDependentModel/eulerianWallDependentModel.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

namespace Foam
{

class eulerianPhasePair;

/*---------------------------------------------------------------------------*\
                      Class eulerianWallDampingModel Declaration
\*---------------------------------------------------------------------------*/

class eulerianWallDampingModel
:
    public eulerianWallDependentModel
{
protected:

    // Protected data

        //- Phase pair
        const eulerianPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("wallDampingModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianWallDampingModel,
            dictionary,
            (
                const dictionary& dict,
                const eulerianPhasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from components
        eulerianWallDampingModel
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~eulerianWallDampingModel();


    // Selectors

        static autoPtr<eulerianWallDampingModel> New
        (
            const word& modelType,
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    // Member Functions

        //- Return damped coefficient
        virtual tmp<volScalarField> damp
        (
            const tmp<volScalarField>&
        ) const = 0;

        //- Return damped force
        virtual tmp<volVectorField> damp
        (
            const tmp<volVectorField>&
        ) const = 0;

        //- Return damped face force
        virtual tmp<surfaceScalarField> damp
        (
            const tmp<surfaceScalarField>&
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

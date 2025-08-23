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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::phaseChangeModel

Description
    Base class for surface film phase change models

SourceFiles
    phaseChangeModel.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "surfaceFilmModels/submodels/filmSubModelBase.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/Fields/scalarField/scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
:
    public filmSubModelBase
{
protected:

    // Protected Member Functions

        //- Latest mass transfer due to phase change
        scalar latestMassPC_;

        //- Total mass transfer due to phase change
        scalar totalMassPC_;


public:

    //- Runtime type information
    TypeName("phaseChange");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        dictionary,
        (
           surfaceFilmRegionModel& film,
           const dictionary& dict
        ),
        (film, dict)
    );


    // Constructors

        //- Construct null
        phaseChangeModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        phaseChangeModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        phaseChangeModel(const phaseChangeModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<phaseChangeModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseChangeModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const scalar dt,
                scalarField& availableMass,
                volScalarField& dMass,
                volScalarField& dEnergy
            );

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            ) = 0;


        // I-O

            //- Provide some feedback
            virtual void info(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2013-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::viscosityModel

Description
    Base class for surface film viscosity models

SourceFiles
    filmViscosityModel.C
    filmViscosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmViscosityModel_H
#define filmViscosityModel_H

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
                     Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
:
    public filmSubModelBase
{
protected:

    // Protected Member Data

        //- Reference to the viscosity field
        volScalarField& mu_;


public:

    //- Runtime type information
    TypeName("viscosity");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        viscosityModel,
        dictionary,
        (
           surfaceFilmRegionModel& film,
           const dictionary& dict,
           volScalarField& mu
        ),
        (film, dict, mu)
    );


    // Constructors

        //- Construct from type name, dictionary and surface film model
        viscosityModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construction
        viscosityModel(const viscosityModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<viscosityModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            volScalarField& mu
        );


    //- Destructor
    virtual ~viscosityModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const volScalarField& p,
                const volScalarField& T
            ) = 0;


        // I-O

            //- Provide some feedback
            virtual void info(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2013-2021 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::turbulenceModel

Description
    Base class for film turbulence models

SourceFiles
    filmTurbulenceModel.C
    filmTurbulenceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmTurbulenceModel_H
#define filmTurbulenceModel_H

#include "surfaceFilmModels/submodels/filmSubModelBase.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class turbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModel
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("turbulenceModel");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        turbulenceModel,
        dictionary,
        (
           surfaceFilmRegionModel& film,
           const dictionary& dict
        ),
        (film, dict)
    );


    // Constructors

        //- Construct null
        turbulenceModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        turbulenceModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        turbulenceModel(const turbulenceModel&) = delete;


    // Selectors

        //- Return a reference to the selected ejection model
        static autoPtr<turbulenceModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulenceModel();


    // Member Functions

        //- Return the film surface velocity
        virtual tmp<volVectorField::Internal> Us() const = 0;

        //- Correct/update the model
        virtual void correct() = 0;

        //- Return the source for the film momentum equation
        virtual tmp<fvVectorMatrix> Su(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

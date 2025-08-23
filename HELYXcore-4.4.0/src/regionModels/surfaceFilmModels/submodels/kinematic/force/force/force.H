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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::force

Description
    Base class for film (stress-based) force models

SourceFiles
    force.C
    forceNew.C

\*---------------------------------------------------------------------------*/

#ifndef force_H
#define force_H

#include "surfaceFilmModels/submodels/filmSubModelBase.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class force Declaration
\*---------------------------------------------------------------------------*/

class force
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("force");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             force,
             dictionary,
             (
                surfaceFilmRegionModel& film,
                const dictionary& dict
             ),
             (film, dict)
         );


    // Constructors

        //- Construct null
        force(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        force
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        force(const force&) = delete;


    // Selectors

        //- Return a reference to the selected force model
        static autoPtr<force> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            const word& modelType
        );


    //- Destructor
    virtual ~force();


    // Member Functions

        //- Correct
        virtual tmp<fvVectorMatrix> correct(volVectorField& U) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const force&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

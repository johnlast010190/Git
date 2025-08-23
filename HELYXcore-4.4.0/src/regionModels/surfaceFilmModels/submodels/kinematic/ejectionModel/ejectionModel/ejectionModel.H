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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::ejectionModel

Description
    Base class for film ejection models, handling mass transfer from the
    film.

SourceFiles
    ejectionModel.C
    ejectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ejectionModel_H
#define ejectionModel_H

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
                      Class ejectionModel Declaration
\*---------------------------------------------------------------------------*/

class ejectionModel
:
    public filmSubModelBase
{
    // Private data

        //- Ejected mass
        scalar ejectedMass_;


protected:

    // Protected Member Functions

        //- Add to ejected mass
        void addToEjectedMass(const scalar dMass);

        //- Correct
        void correct();


public:

    //- Runtime type information
    TypeName("ejectionModel");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ejectionModel,
        dictionary,
        (
           surfaceFilmRegionModel& film,
           const dictionary& dict
        ),
        (film, dict)
    );


    // Constructors

        //- Construct null
        ejectionModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        ejectionModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ejectionModel(const ejectionModel&) = delete;


    // Selectors

        //- Return a reference to the selected ejection model
        static autoPtr<ejectionModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            const word& mdoelType
        );


    //- Destructor
    virtual ~ejectionModel();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToEject,
            scalarField& diameterToEject
        ) = 0;

        //- Return the total mass ejected
        virtual scalar ejectedMassTotal() const;

        //- Accumulate the total mass ejected for the patches into the
        //  scalarField provided
        virtual void patchEjectedMassTotals(scalarField& patchMasses) const
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ejectionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

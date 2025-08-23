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
    Foam::regionModels::surfaceFilmModels::heatTransferModel

Description
    Base class for film heat transfer models

SourceFiles
    heatTransferModel.C
    heatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

#include "surfaceFilmModels/submodels/filmSubModelBase.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("heatTransfer");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        heatTransferModel,
        dictionary,
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        ),
        (film, dict)
    );

    // Constructors

        //- Construct null
        heatTransferModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        heatTransferModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        heatTransferModel(const heatTransferModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<heatTransferModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct() = 0;

            //- Return the heat transfer coefficient [W/m2/K]
            virtual tmp<volScalarField::Internal> h() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heatTransferModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

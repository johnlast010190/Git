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
    Foam::regionModels::surfaceFilmModels::surfaceFilmModel

Description
    Base class for surface film models

SourceFiles
    surfaceFilmModelI.H
    surfaceFilmModel.C
    surfaceFilmModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilmModel_H
#define surfaceFilmModel_H

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFilmModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilmModel
{
public:

    //- Runtime type information
    TypeName("surfaceFilmModel");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceFilmModel,
        mesh,
        (
           const word& modelType,
           const fvMesh& mesh,
           const dimensionedVector& g,
           const word& regionType
        ),
        (modelType, mesh, g, regionType)
    );


    // Constructors

        //- Constructor
        surfaceFilmModel();

        //- Disallow default bitwise copy construction
        surfaceFilmModel(const surfaceFilmModel&) = delete;


    // Selectors

        //- Return a reference to the selected surface film model
        static autoPtr<surfaceFilmModel> New
        (
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType = "surfaceFilm"
        );


    //- Destructor
    virtual ~surfaceFilmModel();


    // Member Functions

        //- Courant number evaluation
        virtual scalar CourantNumber() const = 0;


        // Evolution

            //- Main driver routing to evolve the region - calls other evolves
            virtual void evolve() = 0;


        // Primary region source fields

            //- Return total mass source - Eulerian phase only
            virtual tmp<volScalarField::Internal> Srho() const = 0;

            //- Return mass source for specie i - Eulerian phase only
            virtual tmp<volScalarField::Internal> Srho(const label i) const = 0;

            //- Return momentum source - Eulerian phase only
            virtual tmp<volVectorField::Internal> SU() const = 0;

            //- Return enthalpy source - Eulerian phase only
            virtual tmp<volScalarField::Internal> Sh() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilmModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

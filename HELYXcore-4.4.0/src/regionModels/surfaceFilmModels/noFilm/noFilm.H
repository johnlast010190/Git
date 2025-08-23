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
    Foam::regionModels::surfaceFilmModels::noFilm

Description
    Dummy surface film model for 'none'

SourceFiles
    noFilm.C

\*---------------------------------------------------------------------------*/

#ifndef noFilm_H
#define noFilm_H

#include "surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class noFilm Declaration
\*---------------------------------------------------------------------------*/

class noFilm
:
    public surfaceFilmModel
{
    // Private member data

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noFilm
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType
        );

        //- Disallow default bitwise copy construction
        noFilm(const noFilm&) = delete;


    //- Destructor
    virtual ~noFilm();

    // Member Functions

        //- Courant number evaluation
        virtual scalar CourantNumber() const;

        //- Return total mass source - Eulerian phase only
        virtual tmp<volScalarField::Internal> Srho() const;

        //- Return mass source for specie i - Eulerian phase only
        virtual tmp<volScalarField::Internal> Srho(const label i) const;

        //- Return momentum source - Eulerian phase only
        virtual tmp<volVectorField::Internal> SU() const;

        //- Return enthalpy source - Eulerian phase only
        virtual tmp<volScalarField::Internal> Sh() const;

        //- Main driver routing to evolve the region - calls other evolves
        virtual void evolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

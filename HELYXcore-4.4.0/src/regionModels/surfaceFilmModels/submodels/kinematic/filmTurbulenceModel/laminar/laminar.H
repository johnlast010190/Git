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
    Foam::regionModels::surfaceFilmModels::laminar

Description
    Film laminar turbulence model.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef regionModels_surfaceFilmModels_laminar_H
#define regionModels_surfaceFilmModels_laminar_H

#include "surfaceFilmModels/submodels/kinematic/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                           Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public turbulenceModel
{
private:

    // Private Data

        //- Surface roughness coefficient
        scalar Cf_;


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from surface film model
        laminar(surfaceFilmRegionModel& film, const dictionary& dict);

        //- Disallow default bitwise copy construction
        laminar(const laminar&) = delete;


    //- Destructor
    virtual ~laminar();


    // Member Functions

        //- Return the film surface velocity
        virtual tmp<volVectorField::Internal> Us() const;

        //- Correct/update the model
        virtual void correct();

        //- Return the source for the film momentum equation
        virtual tmp<fvVectorMatrix> Su(volVectorField& U) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laminar&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

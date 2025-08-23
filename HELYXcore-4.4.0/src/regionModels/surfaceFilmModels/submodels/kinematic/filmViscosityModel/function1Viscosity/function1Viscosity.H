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
    (c) 2017-2020 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::function1Viscosity

Description
    The Foam::Function1 temperature-dependent viscosity model multiplies the
    viscosity of a base model by a Function1 temperature function.

SourceFiles
    function1Viscosity.C

\*---------------------------------------------------------------------------*/

#ifndef function1Viscosity_H
#define function1Viscosity_H

#include "surfaceFilmModels/submodels/kinematic/filmViscosityModel/filmViscosityModel/filmViscosityModel.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class function1Viscosity Declaration
\*---------------------------------------------------------------------------*/

class function1Viscosity
:
    public viscosityModel
{
protected:

    // Protected data

        //- Base viscosity model
        autoPtr<viscosityModel> viscosity_;

        //- Viscosity factor as a function of temperature
        autoPtr<Function1<scalar>> function_;


public:

    //- Runtime type information
    TypeName("function1");


    // Constructors

        //- Construct from surface film model
        function1Viscosity
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construct
        function1Viscosity(const function1Viscosity&) = delete;


    //- Destructor
    virtual ~function1Viscosity();


    // Member Functions

        //- Correct
        virtual void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const function1Viscosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

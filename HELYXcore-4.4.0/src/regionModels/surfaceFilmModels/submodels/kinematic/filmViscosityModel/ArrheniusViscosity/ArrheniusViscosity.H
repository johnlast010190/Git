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
    (c) 2015-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::ArrheniusViscosity

Description
    The Arrhenius temperature-dependent viscosity model multiplies the viscosity
    of a base model by an Arrhenius-type temperature function:

        mu = mu0*exp(k1*(1/(T + k2) - 1/(Tref + k2)))

    Where:
        mu0 is the base-model viscosity
        k1 and k2 are Arrhenius coefficients
        T is the local temperature
        Tref is the reference temperature

SourceFiles
    ArrheniusViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusViscosity_H
#define ArrheniusViscosity_H

#include "surfaceFilmModels/submodels/kinematic/filmViscosityModel/filmViscosityModel/filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class ArrheniusViscosity Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusViscosity
:
    public viscosityModel
{
protected:

    // Protected data

        //- Base viscosity model
        autoPtr<viscosityModel> viscosity_;

        //- Coefficient k1
        dimensionedScalar k1_;

        //- Coefficient k2
        dimensionedScalar k2_;

        //- Reference temperature
        dimensionedScalar Tref_;


public:

    //- Runtime type information
    TypeName("Arrhenius");


    // Constructors

        //- Construct from surface film model
        ArrheniusViscosity
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construction
        ArrheniusViscosity(const ArrheniusViscosity&) = delete;


    //- Destructor
    virtual ~ArrheniusViscosity();


    // Member Functions

        //- Correct
        virtual void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ArrheniusViscosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

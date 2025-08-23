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
    Foam::regionModels::surfaceFilmModels::constantViscosity

Description
    Constant viscosity model

SourceFiles
    constantViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef constantViscosity_H
#define constantViscosity_H

#include "surfaceFilmModels/submodels/kinematic/filmViscosityModel/filmViscosityModel/filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class constantViscosity Declaration
\*---------------------------------------------------------------------------*/

class constantViscosity
:
    public viscosityModel
{
protected:

    // Protected data

        //- Constant viscosity [Pa.s]
        dimensionedScalar mu0_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from surface film model
        constantViscosity
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construction
        constantViscosity(const constantViscosity&) = delete;


    //- Destructor
    virtual ~constantViscosity();


    // Member Functions

        //- Correct
        virtual void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const constantViscosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
